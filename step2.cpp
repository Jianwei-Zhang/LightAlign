#include"f.h"
#include"globals.h"
#include<iostream>
#include<fstream>
#include<algorithm>
#include<array>
#include<string>
#include<filesystem> // Generate a folder.
#include <thread>
#include <cstdint>
#include <vector>
#include <cctype>
#include <unordered_map>
#include <cstdio>

using namespace std;

string reverseComplement(string DNAsequence) 
{
	reverse(DNAsequence.begin(), DNAsequence.end());

	static const array<char, 256> complement_map = []() 
		{
		array<char, 256> map{};
		for (int i = 0; i < 256; ++i) map[i] = i;

		map['A'] = 'T'; map['T'] = 'A';
		map['C'] = 'G'; map['G'] = 'C';
		map['a'] = 'T'; map['t'] = 'A';
		map['c'] = 'G'; map['g'] = 'C';
		return map;
		}();

	for (char& b : DNAsequence) 
	{
		b = complement_map[static_cast<unsigned char>(b)];
	}
	return DNAsequence;
}


tuple<int, int, int> MAIN_step2(string dir_name, string inputfile_step2, int groupSize_step3_step4)
{
	int lastGroupSize_step3_step4;
	int groups_step3_step4;
	int l = 0;

	vector<string> block_ID;
	vector<string> block;
	string line;
	string sequenceName;


	ifstream inputFile_0;
	inputFile_0.open(inputfile_step2, ios::binary);

	if (!inputFile_0.is_open())
	{
		std::cerr << "Failed to open " << inputfile_step2 << std::endl;
	}
	size_t pos;
	int FASTA_lines = 0;
	
	while (getline(inputFile_0, line))
	{
		FASTA_lines++;
	}
	
	cout << "total: " << FASTA_lines << endl;
	inputFile_0.close();

	//calculating how many sets are needed
	groups_step3_step4 = (FASTA_lines/2) / groupSize_step3_step4 + 1;
	cout << "groups: " << groups_step3_step4 << endl;
	lastGroupSize_step3_step4 = (FASTA_lines/2) % groupSize_step3_step4;
	cout << "lastGroupsize: " << lastGroupSize_step3_step4 << endl;

	// Create a folder for storing grouped reads
	string grouped_reads_folder_name = dir_name + "/grouped_reads";
	if (std::filesystem::create_directory(grouped_reads_folder_name))
	{
		std::cout << "Successfully created " << grouped_reads_folder_name << std::endl;
	}
	else
	{
		std::cerr << grouped_reads_folder_name << "Creation failed (possibly already exists)" << std::endl;
	}


	int i = 0; // Current group number
	int u = 0; // Line number in group
	int SIZE; // size of current group


	vector<string> ReadsAndID;
	string grouped_reads_file_name;
	SIZE = groupSize_step3_step4;
	ReadsAndID.resize(2 * SIZE);


	grouped_reads_file_name = grouped_reads_folder_name + "/" + to_string(i) + ".txt";
	ofstream outfile1;
	outfile1.open(grouped_reads_file_name);


	ifstream inputFile;
	inputFile.open(inputfile_step2);
	if (!inputFile.is_open()) 
	{
		std::cerr << inputfile_step2 << "Failed to open!" << std::endl;
	}
	while (getline(inputFile, line))
	{
		pos = line.find_first_of("\r\n");
		if (pos != std::string::npos)
		{
			line.resize(pos);
		}
		ReadsAndID[u] = line;

		u++;
		if (u == SIZE * 2)
		{
			//cout << "groups:" << i << endl;
			for (string x : ReadsAndID)
				outfile1 << x << endl;
			outfile1.close();

			i++;
			if (i == groups_step3_step4 - 1)
			{
				SIZE = lastGroupSize_step3_step4;
				ReadsAndID.resize(2 * SIZE);
			}
			else if (i < groups_step3_step4 - 1)
			{
				SIZE = groupSize_step3_step4;
				ReadsAndID.resize(2 * SIZE);
			}
			else
				break;
			u = 0;
			if (i < groups_step3_step4)
			{
				grouped_reads_file_name = grouped_reads_folder_name + "/" + to_string(i) + ".txt"; // Open the file for the next group
				outfile1.open(grouped_reads_file_name);
			}
		}
	}
	inputFile.close();

	
	vector<string> Reads;
	vector<string> ReadsID;
	Reads.resize(SIZE);
	ReadsID.resize(SIZE);

	ofstream outfile(inputfile_step2, std::ios::app);

	//Open the file and start reading.
	for (int i = groups_step3_step4 - 1; i > -1; i--)
	{
		//cout << "groups: " << i << endl;
		grouped_reads_file_name = grouped_reads_folder_name + "/" + to_string(i) + ".txt";
		inputFile.open(grouped_reads_file_name);
		grouped_reads_file_name = grouped_reads_folder_name + "/" + to_string(2* groups_step3_step4 - 1 - i) + ".txt";
		outfile1.open(grouped_reads_file_name, ios::binary);

		if (i == groups_step3_step4 - 1)
		{
			block.resize(lastGroupSize_step3_step4);
			block_ID.resize(lastGroupSize_step3_step4);
		}
		else
		{
			block.resize(groupSize_step3_step4);
			block_ID.resize(groupSize_step3_step4);
		}

		l = 0;
		while (getline(inputFile, line))
		{
			if (l % 2 == 1)
				block[l / 2] = line;
			else
				block_ID[l / 2] = line;

			l = l + 1;
		}
		inputFile.close();

		reverse(block.begin(), block.end());
		reverse(block_ID.begin(), block_ID.end());

		l = 0;
		for (string x : block)
		{
			outfile1 << block_ID[l] << " --complementary" << endl;
			outfile << block_ID[l] << " --complementary" << endl;
			outfile1 << reverseComplement(x) << endl;
			outfile << reverseComplement(x) << endl;
			l++;
		}
		outfile1.close();
	}
	outfile.close();

	// Line offset
	ofstream index_file(dir_name + "/index.txt");
	index_file << "0" << endl;
	inputFile_0.open(inputfile_step2, ios::binary);

	if (!inputFile_0.is_open())
	{
		std::cerr << inputfile_step2 << "Fialed to open!" << std::endl;
	}

	while (getline(inputFile_0, line))
	{
		index_file << static_cast<uint64_t>(inputFile_0.tellg()) << endl;
	}

	inputFile_0.close();
	index_file.close();

	return { FASTA_lines, lastGroupSize_step3_step4, groups_step3_step4 };
}