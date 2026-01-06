#include"f.h"
#include<string>
#include<vector>
#include<iostream>
#include<fstream>
#include<mutex>
#include<unordered_map> // creating hash table
#include<sstream> // Split the string based on commas
#include<tuple>
#include<cmath>
#include <limits>
#include <ctime>
#include<algorithm>
#include <cstdint>
#include <bitset>
#include <cctype>
#include<set>
#include "MurmurHash3.h"
#include"globals.h"
#include <omp.h>
#include<cstdlib>
//#include <intrin.h>

mutex fileMutex;  // Global mutex lock
using namespace std;
const int WORD_SIZE = 64; //Adjust according to machine word length (typically 64-bit)

int8_t match = 1; // 1 byte
int8_t mismatch_penalty = -2;
int8_t gap = -1;

string trim(const string& str) {
	auto start = str.find_first_not_of(" \t\r\n");
	auto end = str.find_last_not_of(" \t\r\n");
	if (start == string::npos) return "";
	return str.substr(start, end - start + 1);
}

string reverse_cigar(const string& cigar) 
{
	vector<string> parts;
	string current_num;

	for (char c : cigar) 
	{
		if (isdigit(c)) 
		{
			current_num += c;
		}
		else 
		{
			parts.push_back(current_num + c);
			current_num.clear();
		}
	}

	std::reverse(parts.begin(), parts.end());

	string reversed;
	for (const string& part : parts) 
	{
		reversed += part;
	}

	return reversed;
}


tuple<string, string, int> needlemanWunsch(const string& seqA, const string& seqB,
	int8_t match, int8_t mismatch, int8_t gap) // seqA/alignA: firstRead; seqB/alignB: secondRead
{
	int actualOverlap = seqB.length();
	int threshold = (1 - error_rate) * match * actualOverlap + error_rate * mismatch * actualOverlap;
	double error_rate_s = error_rate;
	const int segment_size = 1100;

	if (actualOverlap > segment_size)
	{
		error_rate_s = error_rate * 2;

		string alignA_total, alignB_total;
		int total_score = 0;
		int processed = 0;

		alignA_total.reserve(actualOverlap * 1.3);
		alignB_total.reserve(actualOverlap * 1.3);

		while (processed < actualOverlap)
		{
			int seg_len = min(segment_size, actualOverlap - processed);

			string segA = seqA.substr(processed, seg_len);
			string segB = seqB.substr(processed, seg_len);

			auto [alignA_seg, alignB_seg, score_seg] = needlemanWunsch(segA, segB, match, mismatch, gap);

			processed += seg_len;
			total_score += score_seg;

			if (total_score + (actualOverlap - processed) * match < threshold)
			{
				return { "", "", -10000 };
			}

			alignA_total += alignA_seg;
			alignB_total += alignB_seg;
		}

		return { alignA_total, alignB_total, total_score };
	}

	const int m = actualOverlap;
	const int n = actualOverlap;
	const int band_width = max(1, static_cast<int>(m * error_rate_s));

	vector<vector<int>> score_band(m + 1);
	vector<vector<char>> trace_band(m + 1);

	int start_j;
	int end_j;
	int cols;
	int k;
	for (int i = 0; i < m + 1; ++i)
	{
		start_j = max(0, static_cast<int>(i) - band_width);
		end_j = min(static_cast<int>(n), static_cast<int>(i) + band_width);
		cols = end_j - start_j + 1;

		score_band[i].resize(cols, numeric_limits<int>::min());
		trace_band[i].resize(cols, ' ');

		if (i == 0)
		{
			for (int j = start_j; j < 1 + end_j; ++j)
			{
				k = j - start_j;
				score_band[0][k] = j * gap;
			}
		}
	}

	for (size_t i = 1; i < 1 + m; ++i)
	{
		const int start_j = max(0, static_cast<int>(i) - band_width);
		if (start_j == 0)
		{ 
			score_band[i][0] = i * gap;
			trace_band[i][0] = 'U';
		}
	}

	
	int max_score;

	for (int i = 1; i < 1 + m; ++i)
	{
		const int start_j = max(1, static_cast<int>(i) - band_width); 
		const int end_j = min(static_cast<int>(n), static_cast<int>(i) + band_width);
		const int prev_start = max(0, static_cast<int>(i - 1) - band_width); 

		for (int j = start_j; j <= end_j; ++j) 
		{
			max_score = numeric_limits<int>::min();
			char direction = ' ';

			const int curr_start = max(0, static_cast<int>(i) - band_width);
			const int k = j - curr_start;

			const int prev_j = j - 1;
			if (prev_j >= prev_start && prev_j <= min(static_cast<int>(n), static_cast<int>(i - 1) + band_width))
			{
				const int prev_k = prev_j - prev_start;
				const int diag = score_band[i - 1][prev_k] + (seqA[i - 1] == seqB[j - 1] ? match : mismatch);
				if (diag > max_score)
				{
					max_score = diag;
					direction = 'D';
				}
			}

			if (j >= (static_cast<int>(i - 1) - band_width) && j <= (static_cast<int>(i - 1) + band_width))
			{
				const int prev_k = j - prev_start;
				if (prev_k >= 0 && prev_k < static_cast<int>(score_band[i - 1].size()))
				{
					const int up = score_band[i - 1][prev_k] + gap;
					if (up > max_score)
					{
						max_score = up;
						direction = 'U';
					}
				}
			}

			if (k > 0)
			{ 
				const int left = score_band[i][k - 1] + gap;
				if (left > max_score)
				{
					max_score = left;
					direction = 'L';
				}
			}

			score_band[i][k] = max_score;
			trace_band[i][k] = direction;
		}
	}

	// Backtrack to find the optimal path
	string alignA = "";
	string alignB = "";
	int i = m;
	int j = n;
	if (cigar == true)
	{
		while (i > 0 || j > 0)
		{
			const int curr_start = max(0, static_cast<int>(i) - band_width);
			const int k = j - curr_start;

			if (i == 0)
			{
				alignA += '-';
				alignB += seqB[--j];
			}
			else if (j == 0)
			{
				alignA += seqA[--i];
				alignB += '-';
			}
			else
			{
				const char dir = trace_band[i][k];
				switch (dir)
				{
					case 'D':
						alignA += seqA[--i];
						alignB += seqB[--j];
						break;
					case 'U':
						alignA += seqA[--i];
						alignB += '-';
						break;
					case 'L':
						alignA += '-';
						alignB += seqB[--j];
						break;
					default:
						if (i > 0 && j > 0)
						{
							alignA += seqA[--i];
							alignB += seqB[--j];
						}
						else if (i > 0)
						{
							alignA += seqA[--i];
							alignB += '-';
						}
						else
						{
							alignA += '-';
							alignB += seqB[--j];
						}
						break;
				}
		}
	}

	reverse(alignA.begin(), alignA.end());
	reverse(alignB.begin(), alignB.end());
	}
	
	const int final_k = n - (max(0, static_cast<int>(m) - band_width));
	const int final_score = score_band[m][final_k];

	return { alignA, alignB, final_score };
}


bool contains(const std::string& mainStr, const std::string& subStr) // to see if "--complementary" is in read_ID_line
{
	return mainStr.find(subStr) != std::string::npos;
}


// typedef unsigned long int       uint64_t; typedef unsigned long long int  uint64_t;
uint64_t generatingKey(const std::string& AligningUnit)
{
	uint64_t hash[2]; // 128-bit output
	MurmurHash3_x64_128(AligningUnit.data(), AligningUnit.size(), 0, hash);
	return hash[0]; // Use the first 64 bits as the key value
}

tuple<int, int, int, string> difference(const string& str1, const string& str2)  
{
	int mismatch_count = 0; 
	int gap_count = 0; 

	const int total = static_cast<int>(str1.size());
	if (total <= 0) return { 0, 0, 0, "" };

	ostringstream cigar;
	char current_op = ' ';
	int count = 0;
	int m_total = 0;

	auto classify = [](char c1, char c2) 
	{
		if (c1 != '-' && c2 != '-') return 'M';
		return c1 == '-' ? 'D' : 'I';
	};

	for (int i = 0; i < total; ++i) 
	{
		const char op = classify(str1[i], str2[i]);

		if (op == 'M' && str1[i] != str2[i])
			mismatch_count++;
		else if (op == 'I' || op == 'D')
			gap_count++;

		if (op != current_op) 
		{
			if (count > 0) 
			{
				cigar << count << current_op;
			}
			current_op = op;
			count = 1;
		}
		else 
			++count;

		if (op == 'M') ++m_total;
	}

	if (count > 0) cigar << count << current_op;

	return { mismatch_count + gap_count, m_total, total, cigar.str() };
}


string extractID(const string& str) 
{
	string str2 = str.substr(1); //remove @ or >
	size_t space_pos = str2.find(' ');
	if (space_pos != std::string::npos)
	{
		return str2.substr(0, space_pos);
	}
	return str2;
}


void calculateBytesOfBlocks(const std::string& file,
	vector<long long>& groupSizes_bytes,
	int groupSize,
	int groups,
	int lastGroupSize)
{
	int loc = 0;
	int group = 0;
	long long size_sum = 0;
	size_t lineByteSize;
	string line;

	ifstream inputFile4(file, std::ios::binary);

	if (!inputFile4.is_open()) {
		std::cerr << "Failed to open the file: " << file << std::endl;
	}

	while (getline(inputFile4, line)) 
	{
		// Count the number of bytes per line: including inline characters and newlines
		lineByteSize = line.size() + 1;
		size_sum += lineByteSize;
		loc++;

		if (group == groups - 1 || group == groups) 
		{
			if (loc == 2 * lastGroupSize) 
			{
				groupSizes_bytes.push_back(size_sum);
				group++;
				loc = 0;
			}
		}
		else 
		{
			if (loc == 2 * groupSize) 
			{
				groupSizes_bytes.push_back(size_sum);
				group++;
				loc = 0; 
			}
		}
	}

	inputFile4.close();
	std::cout << groupSizes_bytes.size() << endl;
};


//Input a string storing AlignmentUnit 1, and output the splited AlignmentUnit 1.
void seperate(
	vector<uint64_t>& row,
	vector <string>& blockxy_AligningUnit1AndItsLocation,
	vector<vector<uint64_t>>& keyxy_evenAligningUnit,
	vector<string>& blockxy_location_of_evenAligningUnit,
	int groupSize)
{
	int order_in_a_read;
	string buffer;
	uint64_t key;
	string l;
	int a_or_l = 0; // When reading as a AlignmentUnit, determine whether it is an aligningUnit or a location
	
	istringstream one_Read;
	int evenAligningUnit_order;
	int location_of_evenAligningUnit_order;

	//read "blockxy_AligningUnit1AndItsLocation" into "blockxy_evenAligningUnit", "blockxy_location_of_evenAligningUnit"
	for (int i = 0; i < groupSize * 2; i++)
	{
		l = blockxy_AligningUnit1AndItsLocation[i];
	
		a_or_l = i % 2; // aligningUnit or location
		if (a_or_l == 0)
		{
			evenAligningUnit_order = i / 2;

			one_Read.clear(); // Clear status flag
			one_Read.str(l);

			order_in_a_read = 0;
			while (getline(one_Read, buffer, ','))// && order_in_a_read < howManyAlignmentUnits)
			{
				key = stoull(buffer);
				row.push_back(key);
				order_in_a_read++;
			}
			keyxy_evenAligningUnit[evenAligningUnit_order] = row; // The entire row is added to the two-dimensional vector "evenAligningUnit"
			row.clear();
		}

		else
		{
			location_of_evenAligningUnit_order = i / 2;
			blockxy_location_of_evenAligningUnit[location_of_evenAligningUnit_order] = l;
		}
	}
}


// The input stores the string of AlignmentUnit2, and the output is segmented AlignmentUnit2
void seperate2(
	vector<uint64_t>& row,
	vector <string>& blockxy_AligningUnit2AndItsLocation,
	vector<vector<uint64_t>>& keyxy_consecutiveAligningUnit,
	vector<string>& blockxy_location_of_consecutiveAligningUnit,
	int groupSize)
{
	istringstream is;
	string buffer;
	uint64_t key;
	string l;
	for (int i = 0; i < groupSize * 2; i++)
	{
		string l = blockxy_AligningUnit2AndItsLocation[i];
		int u = i % 2;

		if (u == 0)
		{
			int consecutiveAligningUnit_order = i / 2;

			is.clear();
			is.str(l);

			int order_in_a_read = 0;
			while (getline(is, buffer, ',') && order_in_a_read < 2 * DBA_2)
			{
				key = stoull(buffer);
				row[order_in_a_read] = key;
				order_in_a_read++;
			}
			keyxy_consecutiveAligningUnit[consecutiveAligningUnit_order] = row;
		}

		else
		{
			int location_of_consecutiveAligningUnit_order = i / 2;
			blockxy_location_of_consecutiveAligningUnit[location_of_consecutiveAligningUnit_order] = l;
		}
	}
}


void calculate_partitions(
	int group_size,
	int num_threads,
	std::vector<int>& start_inside1,
	std::vector<int>& end_inside1
)
{
	start_inside1.resize(num_threads);
	end_inside1.resize(num_threads);

	const double area = (group_size * group_size) / (2.0 * num_threads);
	double current_position = 0.0;
	double current_length = 2.0 * group_size;

	end_inside1.back() = group_size;

	for (int i = 0; i < num_threads; ++i)
	{
		start_inside1[i] = static_cast<int>(std::round(current_position));

		const double half_length = current_length / 2.0;
		const double radicand = half_length * half_length - area;
		if (radicand < 0)
		{
			throw std::runtime_error("Invalid parameters cause negative sqrt");
		}

		const double delta = half_length - std::sqrt(radicand);
		current_position += delta;
		current_length = 2.0 * std::sqrt(radicand);

		if (i < num_threads - 1)
		{
			end_inside1[i] = static_cast<int>(std::round(current_position));
		}
	}
}

void overlap_between(
	const vector<string>& block11_location_of_evenAligningUnit,
	const vector<string>& block21_location_of_consecutiveAligningUnit,
	const vector<string>& block22_location_of_consecutiveAligningUnit,
	const vector<vector<uint64_t>>& secondRead_key_2,
	const vector<vector<uint64_t>>& secondRead_key_3,
	unordered_map<uint64_t, vector<vector<string> > >& hash_table,
	int howManyReads,
	int firstGroup, int secondGroup, int groups_step3_step4,
	int start,
	int end,
	int windowNumberOfAlignmentUnit,
	int groupSize,
	int lastGroupSize,
	string output,
	bool is_lastGroup
	)
{
	int serial_number_of_firstRead;
	uint64_t key;

	string aligningUnit_read1;
	string aligningUnit_read2;
	int location1;
	int location2;
	string location;
	string overlap_read1;
	string overlap_read2;

	string firstRead_ID;
	string secondRead_ID;

	//int aligningUnit_order_firstRead;
	vector<uint64_t> key2;
	vector<uint64_t> key3;

	set<string> firstRead_used; // if (firstGroup == secondGroup)
	set<vector<int>>overlapped_set; // if (firstGroup != secondGroup)
	vector<int>overlapped(2);

	ofstream outfile_filted(output, std::ios::app);
	for (int serial_number_of_secondRead = start; serial_number_of_secondRead < end; serial_number_of_secondRead++)
	{
		secondRead_ID = block21_location_of_consecutiveAligningUnit[serial_number_of_secondRead];
		key2 = secondRead_key_2[serial_number_of_secondRead];

		for (int n1 = 0; n1 < 2 * DBA_2; n1++) // // The first alignment units of firstRead and secondRead overlap
		{
			key = key2[n1];
			auto it = hash_table.find(key);
			
			// Pairwise comparison
			if (it != hash_table.end())
			{
				for(vector<string> location : it->second)
				{ 
					serial_number_of_firstRead = stoi(location[0]); // location format：to_string(x), to_string(aligningUnit_order_firstRead * DBA)
					firstRead_ID = block11_location_of_evenAligningUnit[serial_number_of_firstRead];
					location1 = stoi(location[1]); // evenAligningUnit's starting position on read
					
					if (n1 % 2 == 0)
						location2 = n1/2;
					else if (n1 % 2 == 1)
						location2 = -1 * (lengthOfAlignmentUnit + n1/2);

					if (firstRead_ID != secondRead_ID)
					{
						if (firstGroup == secondGroup)
						{
							if (firstRead_used.find(firstRead_ID) != firstRead_used.end())
							{
								continue;
							}
							else
							{
								//outfile_filted << firstRead_ID << "\t" << location1 << "\t" << secondRead_ID << "\t" << location2 << endl;
								overlapped[0] = serial_number_of_firstRead;
								overlapped[1] = serial_number_of_secondRead;
								if (overlapped[0] > overlapped[1]) { swap(overlapped[0], overlapped[1]); }
								if (overlapped_set.find(overlapped) != overlapped_set.end()) // The overlapped read pairs already exists
								{
									continue;
								}
								else
								{
									overlapped_set.insert(overlapped);
									outfile_filted << firstRead_ID << "\t" << location1 << "\t" << secondRead_ID << "\t" << location2 << endl;
								}
							}
						}
						else 
						{
							overlapped[0] = serial_number_of_firstRead;
							overlapped[1] = serial_number_of_secondRead;
							
							if (overlapped[0] > overlapped[1]) {swap(overlapped[0], overlapped[1]);}
							if (overlapped_set.find(overlapped) != overlapped_set.end())
							{
								continue;
							}
							else 
							{
								overlapped_set.insert(overlapped);
								outfile_filted << firstRead_ID << "\t" << location1 << "\t" << secondRead_ID << "\t" << location2 << endl;
							}
						}
					}
				}

			}
		}
		if (firstGroup == secondGroup)
			firstRead_used.insert(secondRead_ID);
	}

	firstRead_used.clear();
	overlapped_set.clear();

	for (int serial_number_of_secondRead = start; serial_number_of_secondRead < end; serial_number_of_secondRead++)
	{
		secondRead_ID = block22_location_of_consecutiveAligningUnit[serial_number_of_secondRead]; // [0] [0] ;
		key3 = secondRead_key_3[serial_number_of_secondRead];
		for (int n1 = 0; n1 < 2 * DBA_2; n1++)
		{
			key = key3[n1];
			auto it = hash_table.find(key);
				
			// Pairwise comparison0
			if (it != hash_table.end())
			{
				for (vector<string> location : it->second)
				{
					serial_number_of_firstRead = stoi(location[0]);
					firstRead_ID = block11_location_of_evenAligningUnit[serial_number_of_firstRead];
					location1 = stoi(location[1]); // evenAligningUnit's starting position on read
					
					if (n1 % 2 == 0)
						location2 = n1/2;
					else if (n1 % 2 == 1)
						location2 = -1 * (lengthOfAlignmentUnit + n1/2);
					
					if (firstRead_ID != secondRead_ID)
					{
						if (firstGroup == secondGroup)
						{
							if (firstRead_used.find(firstRead_ID) != firstRead_used.end()) 
							{
								continue;
							}
							else 
							{
								//outfile_filted << firstRead_ID << "\t" << location1 << "\t" << secondRead_ID << "\t" << location2 << endl;
								overlapped[0] = serial_number_of_firstRead;
								overlapped[1] = serial_number_of_secondRead;
								
								if (overlapped[0] > overlapped[1]) { swap(overlapped[0], overlapped[1]); }
								if (overlapped_set.find(overlapped) != overlapped_set.end())
								{
									continue;
								}
								else
								{
									overlapped_set.insert(overlapped);
									outfile_filted << firstRead_ID << "\t" << location1 << "\t" << secondRead_ID << "\t" << location2 << endl;
								}
							}
						}
						else 
						{
							overlapped[0] = serial_number_of_firstRead;
							overlapped[1] = serial_number_of_secondRead;
							
							if (overlapped[0] > overlapped[1]) { swap(overlapped[0], overlapped[1]); }
							if (overlapped_set.find(overlapped) != overlapped_set.end())
							{
								continue;
							}
							else
							{
								overlapped_set.insert(overlapped);
								outfile_filted << firstRead_ID << "\t" << location1 << "\t" << secondRead_ID << "\t" << location2 << endl;
							}
						}
					}
				}
			}
		}
		if (firstGroup == secondGroup)
			firstRead_used.insert(extractID(secondRead_ID));
	}

}


long long countLinesFast(const string& filename)
{
	ifstream file(filename, ios::binary);
	if (!file.is_open())
	{
		cerr << "Error opening file: " << filename << endl;
		return -1;
	}

	long long lineCount = 0;
	char c;

	while (file.get(c))
	{
		if (c == '\r')
		{
			// Handle cases of \r or \r\n (different line ending formats)
			if (file.peek() == '\n')
			{
				file.get();
			}
			lineCount++;
			if (lineCount % 1000000 == 0)
				cout << lineCount << endl;
		}
		else if (c == '\n')
		{
			lineCount++;
			if (lineCount % 1000000 == 0)
				cout << lineCount << endl;
		}
	}

	// Handle the case where the last line lacks a terminating newline character.
	file.clear();
	file.seekg(-1, ios::end);
	if (file.get(c) && c != '\n' && c != '\r')
	{
		lineCount++;
	}

	file.close();
	return lineCount;
}


// PAF is 1-based
void alignment(string filtedFile,
	ofstream& outfile_overlap,
	ofstream& outfile_score,
	unordered_map<string, long long>& reads_hashtable,
	string outputfile_step1)
{
	int actualOverlap;
	int read1_len;
	int read2_len;
	string read1_id;
	string read2_id;
	string read1;
	string read2;
	string overlap_read1;
	string overlap_read2;
	string firstRead_ID;
	string secondRead_ID;

	int read1_start; //start of overlap region
	int read1_end;
	int read2_start;
	int read2_end;
	int location1;
	int location2;

	size_t space_pos;
	// forward/reverse
	string Plus_or_minus;

	ifstream inputFile1; // filtedFASTA.txt
	inputFile1.open(outputfile_step1, ios::binary);

	ifstream inputFile; // filted.txt
	inputFile.open(filtedFile, ios::binary);

	string line;
	string token;
	vector<string> tokens;
	while (getline(inputFile, line))
	{
		istringstream iss(line);

		while (getline(iss, token, '\t')) 
		{  // splited by TAB
			tokens.push_back(token);
		}
		
		location1 = stoi(tokens[1]);
		location2 = stoi(tokens[3]);
		read1_id = tokens[0];
		read2_id = tokens[2];
		// Find the position of the first space
		firstRead_ID = extractID(read1_id); // no "--complementary"
		secondRead_ID = extractID(read2_id); // no "--complementary"
		if (firstRead_ID == secondRead_ID)
		{
			tokens.clear();
			continue;
		}

		tokens.clear();

		auto it = reads_hashtable.find(read1_id);
		inputFile1.clear();
		inputFile1.seekg(it->second);
		getline(inputFile1, read1);

		it = reads_hashtable.find(read2_id);
		inputFile1.clear();
		inputFile1.seekg(it->second);
		getline(inputFile1, read2);

		read1_len = read1.length();
		read2_len = read2.length();

		if (location2 < 0) // The comparison is to the right of read2
		{
			actualOverlap = location1 - location2;

			if (read1_len > actualOverlap && read2_len > actualOverlap)
			{
				overlap_read1 = read1.substr(0, actualOverlap);
				overlap_read2 = read2.substr(read2_len - actualOverlap, actualOverlap);
				auto [alignmentA, alignmentB, score] = needlemanWunsch(overlap_read1, overlap_read2, match, mismatch_penalty, gap);
				outfile_score << score << "\t" << actualOverlap << endl;

				if (score >= actualOverlap * match * (1 - error_rate) + actualOverlap * mismatch_penalty * error_rate)
				{
					auto [NM, m_total, total, CIGAR] = difference(alignmentA, alignmentB);

					if (contains(read1_id, "--complementary"))
					{
						Plus_or_minus = "-";
						read2_start = read2_len - actualOverlap;
						read2_end = read2_len;

						read1_start = read1_len - actualOverlap;
						read1_end = read1_len;
					}
					else if (contains(read2_id, "--complementary"))
					{
						CIGAR = reverse_cigar(CIGAR);
						Plus_or_minus = "-";
						read1_start = 0;
						read1_end = actualOverlap;

						read2_start = 0;
						read2_end = actualOverlap;
					}
					else
					{
						Plus_or_minus = "+";
						read1_start = 0;
						read1_end = actualOverlap;

						read2_start = read2_len - actualOverlap;
						read2_end = read2_len;
					}
					lock_guard<mutex> lock(fileMutex);
					outfile_overlap << firstRead_ID << "\t"
						<< read1_len << "\t"
						<< read1_start + 1 << "\t"
						<< read1_end << "\t"
						<< Plus_or_minus << "\t"
						<< secondRead_ID << "\t" //location_of_evenAligningUnit[serial_number_of_firstRead][0][0]
						<< read2_len << "\t"
						<< read2_start + 1 << "\t"
						<< read2_end << "\t"
						<< m_total << "\t"
						<< total << "\t"
						<< "0" << "\t"
						<< "cg:Z:" << CIGAR << "\t"
						<< "NM:i:" << NM // Total number of mismatches and gaps in the alignment
						<< endl; // secondRead front
				}
			}
		
			else if (read1_len <= actualOverlap && read2_len > actualOverlap)
			{
				overlap_read1 = read1;
				overlap_read2 = read2.substr(read2_len - actualOverlap, read1_len);
				auto [alignmentA, alignmentB, score] = needlemanWunsch(overlap_read1, overlap_read2, match, mismatch_penalty, gap);
				outfile_score << score << "\t" << actualOverlap << endl;

				if (score >= actualOverlap * match * (1 - error_rate) + actualOverlap * mismatch_penalty * error_rate)
				{
					auto [NM, m_total, total, CIGAR] = difference(alignmentA, alignmentB);

					if (contains(read1_id, "--complementary"))
					{
						Plus_or_minus = "-";
						read2_start = read2_len - actualOverlap;
						read2_end = read2_len - actualOverlap + read1_len;

						read1_start = 0;
						read1_end = read1_len;
					}
					else if (contains(read2_id, "--complementary"))
					{
						CIGAR = reverse_cigar(CIGAR);
						Plus_or_minus = "-";
						read1_start = 0;
						read1_end = read1_len;

						read2_start = actualOverlap - read1_len;
						read2_end = actualOverlap;
					}
					else
					{
						Plus_or_minus = "+";
						read1_start = 0;
						read1_end = read1_len;

						read2_start = read2_len - actualOverlap;
						read2_end = read2_len - actualOverlap + read1_len;
					}
					lock_guard<mutex> lock(fileMutex);
					outfile_overlap << firstRead_ID << "\t"
						<< read1_len << "\t"
						<< read1_start + 1 << "\t"
						<< read1_end << "\t"
						<< Plus_or_minus << "\t"
						<< secondRead_ID << "\t" //location_of_evenAligningUnit[serial_number_of_firstRead][0][0]
						<< read2_len << "\t"
						<< read2_start + 1 << "\t"
						<< read2_end << "\t"
						<< m_total << "\t"
						<< total << "\t"
						<< "0" << "\t"
						<< "cg:Z:" << CIGAR << "\t"
						<< "NM:i:" << NM // Total number of mismatches and gaps in the alignment
						<< endl; // secondRead front
				}
			}

			else if (read1_len > actualOverlap && read2_len <= actualOverlap)
			{
				overlap_read1 = read1.substr(actualOverlap - read2_len, read2_len);
				overlap_read2 = read2;
				auto [alignmentA, alignmentB, score] = needlemanWunsch(overlap_read1, overlap_read2, match, mismatch_penalty, gap);
				outfile_score << score << "\t" << read2_len << endl;

				if (score >= read2_len * match * (1 - error_rate) + read2_len * mismatch_penalty * error_rate)
				{
					auto [NM, m_total, total, CIGAR] = difference(alignmentA, alignmentB);

					if (contains(read1_id, "--complementary"))
					{
						Plus_or_minus = "-";
						read2_start = 0;
						read2_end = read2_len;

						read1_start =  read1_len - actualOverlap;
						read1_end = read1_len - actualOverlap + read2_len;
					}
					else if (contains(read2_id, "--complementary"))
					{
						CIGAR = reverse_cigar(CIGAR);
						Plus_or_minus = "-";
						read1_start = actualOverlap - read2_len;
						read1_end = actualOverlap;

						read2_start = 0;
						read2_end = read2_len;
					}
					else
					{
						Plus_or_minus = "+";
						read1_start = actualOverlap - read2_len;
						read1_end = actualOverlap;

						read2_start = 0;
						read2_end = read2_len;
					}
					lock_guard<mutex> lock(fileMutex);
					outfile_overlap << firstRead_ID << "\t"
						<< read1_len << "\t"
						<< read1_start + 1 << "\t"
						<< read1_end << "\t"
						<< Plus_or_minus << "\t"
						<< secondRead_ID << "\t" //location_of_evenAligningUnit[serial_number_of_firstRead][0][0]
						<< read2_len << "\t"
						<< read2_start + 1 << "\t"
						<< read2_end << "\t"
						<< m_total << "\t"
						<< total << "\t"
						<< "0" << "\t"
						<< "cg:Z:" << CIGAR << "\t"
						<< "NM:i:" << NM // Total number of mismatches and gaps in the alignment
						<< endl; // secondRead front
				}
			}
		
			else if (read1_len <= actualOverlap && read2_len <= actualOverlap)
			{
				actualOverlap = (read2_len + location2) + (read1_len - location1);
				overlap_read1 = read1.substr(read1_len - actualOverlap, actualOverlap);
				overlap_read2 = read2.substr(0, actualOverlap);
				auto [alignmentA, alignmentB, score] = needlemanWunsch(overlap_read1, overlap_read2, match, mismatch_penalty, gap);
				outfile_score << score << "\t" << actualOverlap << endl;

				if (score >= actualOverlap * match * (1 - error_rate) + actualOverlap * mismatch_penalty * error_rate)
				{
					auto [NM, m_total, total, CIGAR] = difference(alignmentA, alignmentB);

					if (contains(read1_id, "--complementary"))
					{
						Plus_or_minus = "-";
						read1_start = 0;
						read1_end = actualOverlap;

						read2_start = 0;
						read2_end = actualOverlap;
					}
					else if (contains(read2_id, "--complementary"))
					{
						CIGAR = reverse_cigar(CIGAR);
						Plus_or_minus = "-";
						read1_start = read1_len - actualOverlap;
						read1_end = read1_len;

						read2_start = read2_len - actualOverlap;
						read2_end = read2_len;
					}
					else
					{
						Plus_or_minus = "+";
						read1_start = read1_len - actualOverlap;
						read1_end = read1_len;

						read2_start = 0;
						read2_end = actualOverlap;
					}

					lock_guard<mutex> lock(fileMutex); // lock before write
					outfile_overlap << firstRead_ID << "\t"
						<< read1_len << "\t"
						<< read1_start + 1 << "\t"
						<< read1_end << "\t"
						<< Plus_or_minus << "\t"
						<< secondRead_ID << "\t"
						<< read2_len << "\t"
						<< read2_start + 1 << "\t"
						<< read2_end << "\t"
						<< m_total << "\t"
						<< total << "\t"
						<< "0" << "\t"
						<< "cg:Z:" << CIGAR << "\t"
						<< "NM:i:" << NM // Total number of mismatches and gaps in the alignment
						<< endl; // secondRead front
				}
			}
		}

		else // The comparison is to the left of read2
		{
			actualOverlap = location2 + (read1_len - location1);

			if (read1_len > actualOverlap && read2_len > actualOverlap)
			{
				overlap_read1 = read1.substr(read1_len - actualOverlap, actualOverlap);
				overlap_read2 = read2.substr(0, actualOverlap);
				auto [alignmentA, alignmentB, score] = needlemanWunsch(overlap_read1, overlap_read2, match, mismatch_penalty, gap);
				outfile_score << score << "\t" << actualOverlap << endl;

				if (score >= actualOverlap * match * (1 - error_rate) + actualOverlap * mismatch_penalty * error_rate)
				{
					auto [NM, m_total, total, CIGAR] = difference(alignmentA, alignmentB);

					if (contains(read1_id, "--complementary"))
					{
						Plus_or_minus = "-";
						read1_start = 0;
						read1_end = actualOverlap;

						read2_start = 0;
						read2_end = actualOverlap;
					}
					else if (contains(read2_id, "--complementary"))
					{
						CIGAR = reverse_cigar(CIGAR);
						Plus_or_minus = "-";
						read1_start = read1_len - actualOverlap;
						read1_end = read1_len;

						read2_start = read2_len - actualOverlap;
						read2_end = read2_len;
					}
					else
					{
						Plus_or_minus = "+";
						read1_start = read1_len - actualOverlap;
						read1_end = read1_len;

						read2_start = 0;
						read2_end = actualOverlap;
					}

					lock_guard<mutex> lock(fileMutex); // lock before write
					outfile_overlap << firstRead_ID << "\t"
						<< read1_len << "\t"
						<< read1_start + 1 << "\t"
						<< read1_end << "\t"
						<< Plus_or_minus << "\t"
						<< secondRead_ID << "\t"
						<< read2_len << "\t"
						<< read2_start + 1 << "\t"
						<< read2_end << "\t"
						<< m_total << "\t"
						<< total << "\t"
						<< "0" << "\t"
						<< "cg:Z:" << CIGAR << "\t"
						<< "NM:i:" << NM // Total number of mismatches and gaps in the alignment
						<< endl; // secondRead front
				}
			}
			
			else if (read1_len <= actualOverlap && read2_len > actualOverlap)
			{
				overlap_read1 = read1;
				overlap_read2 = read2.substr(actualOverlap - read1_len, read1_len);
				auto [alignmentA, alignmentB, score] = needlemanWunsch(overlap_read1, overlap_read2, match, mismatch_penalty, gap);
				outfile_score << score << "\t" << actualOverlap << endl;

				if (score >= actualOverlap * match * (1 - error_rate) + actualOverlap * mismatch_penalty * error_rate)
				{
					auto [NM, m_total, total, CIGAR] = difference(alignmentA, alignmentB);

					if (contains(read1_id, "--complementary"))
					{
						Plus_or_minus = "-";
						read1_start = 0;
						read1_end = read1_len;

						read2_start = actualOverlap - read1_len;
						read2_end = actualOverlap;
					}
					else if (contains(read2_id, "--complementary"))
					{
						CIGAR = reverse_cigar(CIGAR);
						Plus_or_minus = "-";
						read1_start = 0;
						read1_end = read1_len;

						read2_start = read2_len - actualOverlap;
						read2_end = read2_len - actualOverlap + read1_len;
					}
					else
					{
						Plus_or_minus = "+";
						read1_start = 0;
						read1_end = read1_len;

						read2_start = actualOverlap - read1_len;
						read2_end = actualOverlap;
					}

					lock_guard<mutex> lock(fileMutex); // lock before write
					outfile_overlap << firstRead_ID << "\t"
						<< read1_len << "\t"
						<< read1_start + 1 << "\t"
						<< read1_end << "\t"
						<< Plus_or_minus << "\t"
						<< secondRead_ID << "\t"
						<< read2_len << "\t"
						<< read2_start + 1 << "\t"
						<< read2_end << "\t"
						<< m_total << "\t"
						<< total << "\t"
						<< "0" << "\t"
						<< "cg:Z:" << CIGAR << "\t"
						<< "NM:i:" << NM // Total number of mismatches and gaps in the alignment
						<< endl; // secondRead front
				}
			}
			
			else if (read1_len > actualOverlap && read2_len <= actualOverlap)
			{
				overlap_read1 = read1.substr(read1_len - actualOverlap, read2_len);
				overlap_read2 = read2;
				auto [alignmentA, alignmentB, score] = needlemanWunsch(overlap_read1, overlap_read2, match, mismatch_penalty, gap);
				outfile_score << score << "\t" << read2_len << endl;

				if (score >= read2_len * match * (1 - error_rate) + read2_len * mismatch_penalty * error_rate)
				{
					auto [NM, m_total, total, CIGAR] = difference(alignmentA, alignmentB);

					if (contains(read1_id, "--complementary"))
					{
						Plus_or_minus = "-";
						read1_start = actualOverlap - read2_len;
						read1_end = actualOverlap;

						read2_start = 0;
						read2_end = read2_len;
					}
					else if (contains(read2_id, "--complementary"))
					{
						CIGAR = reverse_cigar(CIGAR);
						Plus_or_minus = "-";
						read1_start = read1_len - actualOverlap;
						read1_end = read1_len - actualOverlap + read2_len;

						read2_start = 0;
						read2_end = read2_len;
					}
					else
					{
						Plus_or_minus = "+";
						read1_start = read1_len - actualOverlap;
						read1_end = read1_len - actualOverlap + read2_len;

						read2_start = 0;
						read2_end = read2_len;
					}

					lock_guard<mutex> lock(fileMutex); // lock before write
					outfile_overlap << firstRead_ID << "\t"
						<< read1_len << "\t"
						<< read1_start + 1 << "\t"
						<< read1_end << "\t"
						<< Plus_or_minus << "\t"
						<< secondRead_ID << "\t"
						<< read2_len << "\t"
						<< read2_start + 1 << "\t"
						<< read2_end << "\t"
						<< m_total << "\t"
						<< total << "\t"
						<< "0" << "\t"
						<< "cg:Z:" << CIGAR << "\t"
						<< "NM:i:" << NM // Total number of mismatches and gaps in the alignment
						<< endl; // secondRead front
				}
			}
		
			else if (read1_len <= actualOverlap && read2_len <= actualOverlap) // sort to make sure read1_len > read2_len
			{
				actualOverlap = read2_len - location2;
				overlap_read1 = read1.substr(0, read2_len - location2);
				overlap_read2 = read2.substr(location2, read2_len - location2);
				auto [alignmentA, alignmentB, score] = needlemanWunsch(overlap_read1, overlap_read2, match, mismatch_penalty, gap);
				outfile_score << score << "\t" << actualOverlap << endl;

				if (score >= actualOverlap * match * (1 - error_rate) + actualOverlap * mismatch_penalty * error_rate)
				{
					auto [NM, m_total, total, CIGAR] = difference(alignmentA, alignmentB);

					if (contains(read1_id, "--complementary"))
					{
						Plus_or_minus = "-";
						read1_start = read1_len - actualOverlap;
						read1_end = read1_len;

						read2_start = location2;
						read2_end = read2_len;
					}
					else if (contains(read2_id, "--complementary"))
					{
						CIGAR = reverse_cigar(CIGAR);
						Plus_or_minus = "-";
						read1_start = 0;
						read1_end = actualOverlap;

						read2_start = 0;
						read2_end = actualOverlap;
					}
					else
					{
						Plus_or_minus = "+";
						read1_start = 0;
						read1_end = actualOverlap;

						read2_start = location2;
						read2_end = read2_len;
					}

					lock_guard<mutex> lock(fileMutex); // lock before write
					outfile_overlap << firstRead_ID << "\t"
						<< read1_len << "\t"
						<< read1_start + 1 << "\t"
						<< read1_end << "\t"
						<< Plus_or_minus << "\t"
						<< secondRead_ID << "\t"
						<< read2_len << "\t"
						<< read2_start + 1 << "\t"
						<< read2_end << "\t"
						<< m_total << "\t"
						<< total << "\t"
						<< "0" << "\t"
						<< "cg:Z:" << CIGAR << "\t"
						<< "NM:i:" << NM // Total number of mismatches and gaps in the alignment
						<< endl; // secondRead front
				}
			}
		}
	}
}



// Both of the following functions are used only in firstGroup
void read_in_data_for_allReads(ifstream& inputFile,
	string& dir_name, string& grouped_reads_folder_name,
	vector<string>& block11_allReads, 
	vector<string>& block12_allReads,
	vector<int>& which_group_had_alignment_inside,
	int groupSize, int firstGroup,
	int groupSize_step3_step4, int groups_step3_step4, int lastGroupSize_step3_step4)
{
	string line;
	int loc = 0;
	int loc_inside_11 = 0;
	int loc_inside_12 = 0;

	// Position the file pointer to the "startByte" location
	string filename = dir_name + "/" + grouped_reads_folder_name + "/" + to_string(firstGroup) + ".txt";
	inputFile.open(filename);

	for (int lineCount = 0; getline(inputFile, line) && loc_inside_11 < groupSize; ++lineCount) 
	{
		if (lineCount % 2 == 0) continue;
		block11_allReads[loc_inside_11++] = line;
	}
	inputFile.close();

	if (find(which_group_had_alignment_inside.begin(), which_group_had_alignment_inside.end(), firstGroup) == which_group_had_alignment_inside.end())
	{
		loc = 0;
		filename = dir_name + "/" + grouped_reads_folder_name + "/" + to_string((groups_step3_step4 * 2 - 1) - firstGroup) + ".txt";
		inputFile.open(filename);

		for (int lineCount = 0; getline(inputFile, line) && loc_inside_12 < groupSize; ++lineCount)
		{
			if (lineCount % 2 == 0) continue;
			block12_allReads[loc_inside_12++] = line;
		}
		inputFile.close();
	}
}



void read_in_data_for_all_AligningUnitAndItsLocation(ifstream& inputFile,
	string& dir_name, string& grouped_aligningUnit_folder_name,
	vector<string>& block11_AligningUnitAndItsLocation,
	vector<string>& block12_AligningUnitAndItsLocation,
	vector<int>& which_group_had_alignment_inside,
	int groupSize, int firstGroup,
	int groupSize_step3_step4, int groups_step3_step4, int lastGroupSize_step3_step4) // AligningUnit1 is treated the same as AligningUnit2
{
	string line;
	int loc_inside_11;
	int loc_inside_12;
	loc_inside_11 = 0;
	loc_inside_12 = 0;

	string filename = dir_name + "/" + grouped_aligningUnit_folder_name + "/" + to_string(firstGroup) + ".txt";
	inputFile.open(filename);
	while (getline(inputFile, line))
	{
		if (loc_inside_11 < 2 * groupSize)
		{
			block11_AligningUnitAndItsLocation[loc_inside_11] = line;
			loc_inside_11 = loc_inside_11 + 1;
		}
		else
			break;
	}
	inputFile.close();

	if (find(which_group_had_alignment_inside.begin(), which_group_had_alignment_inside.end(), firstGroup) == which_group_had_alignment_inside.end())
	{
		filename = dir_name + "/" + grouped_aligningUnit_folder_name + "/" + to_string((groups_step3_step4 * 2 - 1) - firstGroup) + ".txt";
		inputFile.open(filename);

		while (getline(inputFile, line))
		{
			if (loc_inside_12 < 2 * groupSize)
			{
				block12_AligningUnitAndItsLocation[loc_inside_12] = line;
				loc_inside_12 = loc_inside_12 + 1;
			}
			else
				break;
		}
		inputFile.close();
	}
}


// Both of the following functions are used only in secondGroup
void read_in_data_for_allReads_secondGroup(ifstream& inputFile,
	string& dir_name, string& grouped_reads_folder_name,
	vector<string>& block21_allReads,
	vector<string>& block22_allReads,
	int groupSize, int secondGroup,
	int groupSize_step3_step4, int groups_step3_step4, int lastGroupSize_step3_step4)
{
	string line;
	int i = 0;
	int loc_inside_21 = 0;
	int loc_inside_22 = 0;

	string filename = dir_name + "/" + grouped_reads_folder_name + "/" + to_string(secondGroup) + ".txt";
	inputFile.open(filename);
	for (int lineCount = 0; getline(inputFile, line) && loc_inside_21 < groupSize; ++lineCount)
	{
		if (lineCount % 2 == 0) continue;
		block21_allReads[loc_inside_21++] = line;
	}
	inputFile.close();

	i = 0;
	filename = dir_name + "/" + grouped_reads_folder_name + "/" + to_string((groups_step3_step4 * 2 - 1) - secondGroup) + ".txt";
	inputFile.open(filename);
	for (int lineCount = 0; getline(inputFile, line) && loc_inside_22 < groupSize; ++lineCount)
	{
		if (lineCount % 2 == 0) continue;
		block22_allReads[loc_inside_22++] = line;
	}
	inputFile.close();
}


void read_in_data_for_all_AligningUnitAndItsLocation_secondGroup(ifstream& inputFile,
	string& dir_name, string& grouped_aligningUnit_folder_name,
	vector<string>& block21_AligningUnitAndItsLocation,
	vector<string>& block22_AligningUnitAndItsLocation,
	int groupSize, int secondGroup, 
	int groupSize_step3_step4, int groups_step3_step4, int lastGroupSize_step3_step4) // AligningUnit1 is treated the same as AligningUnit2
{
	string line;
	int i = 0;
	int loc_inside_21 = 0;
	int loc_inside_22 = 0;

	string filename = dir_name + "/" + grouped_aligningUnit_folder_name + "/" + to_string(secondGroup) + ".txt";
	inputFile.open(filename);
	while (getline(inputFile, line))
	{
		if (loc_inside_21 < 2 * groupSize)
		{
			block21_AligningUnitAndItsLocation[loc_inside_21] = line;
			loc_inside_21 = loc_inside_21 + 1;
		}
		else
			break;
	}
	inputFile.close();

	i = 0;
	filename = dir_name + "/" + grouped_aligningUnit_folder_name + "/" + to_string((groups_step3_step4 * 2 - 1) - secondGroup) + ".txt";
	inputFile.open(filename);
	while (getline(inputFile, line))
	{
		if (loc_inside_22 < 2 * groupSize)
		{
			block22_AligningUnitAndItsLocation[loc_inside_22] = line;
			loc_inside_22 = loc_inside_22 + 1;
		}
		else
			break;
	}
	inputFile.close();
}
