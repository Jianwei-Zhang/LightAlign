#include<iostream>
#include<string>
#include <memory>
#include<fstream>
#include<filesystem> // Generate a folder.
#include <thread>
#include"f.h"
#include"globals.h"

using namespace std;

string get01Version(vector<int>& ATaccumulate, int startpoint, int howManyWindowInAnAlignmentUnit)
{
	
	string fluctuation = "1";
	int former_ATwindow; //AT content of first window
	int ATwindow;
	int translocated_position = startpoint + step - 1;
	int len = ATaccumulate.size();

	if (startpoint == 0)
	{
		former_ATwindow = ATaccumulate[window - 1];
		for (int i = 1; i < howManyWindowInAnAlignmentUnit; i++)
		{
			//translocated_position = 0的时候第一个window其实没算进去，ATaccumulate[0]被减掉了
			ATwindow = ATaccumulate[translocated_position + window] - ATaccumulate[translocated_position];

			if (ATwindow >= former_ATwindow)
				fluctuation = fluctuation + "1";
			else
				fluctuation = fluctuation + "0";

			former_ATwindow = ATwindow;
			translocated_position = translocated_position + step;
		}
	}

	else if (startpoint < 0)
	{
		former_ATwindow = ATaccumulate[len + startpoint + window - 1] - ATaccumulate[len + startpoint - 1];
		for (int i = 1; i < howManyWindowInAnAlignmentUnit; i++)
		{
			//translocated_position = 0的时候第一个window其实没算进去，ATaccumulate[0]被减掉了
			ATwindow = ATaccumulate[len + translocated_position + window] - ATaccumulate[len + translocated_position];

			if (ATwindow >= former_ATwindow)
				fluctuation = fluctuation + "1";
			else
				fluctuation = fluctuation + "0";

			former_ATwindow = ATwindow;
			translocated_position = translocated_position + step;
		}
	}
	
	else // > 0
	{
		former_ATwindow = ATaccumulate[startpoint + window - 1] - ATaccumulate[startpoint - 1];
		for (int i = 1; i < howManyWindowInAnAlignmentUnit; i++)
		{
			//translocated_position = 0的时候第一个window其实没算进去，ATaccumulate[0]被减掉了
			ATwindow = ATaccumulate[translocated_position + window] - ATaccumulate[translocated_position];

			if (ATwindow >= former_ATwindow)
				fluctuation = fluctuation + "1";
			else
				fluctuation = fluctuation + "0";

			former_ATwindow = ATwindow;
			translocated_position = translocated_position + step;
		}
	}

	return fluctuation;
}


void gettingAligningUnitAndItsLocation(
	vector<string>& Reads,
	vector<string>& ReadsID,
	int start,
	int end,
	vector<vector<uint64_t>>& v_e, vector<string>&v_e_ID,
	vector<vector<uint64_t>>& v_c, vector<string>& v_c_ID)
{
	int howManyWindowInAnAlignmentUnit = (lengthOfAlignmentUnit - window) / step + 1;

	int howManyAnAlignmentUnitsInARead;

	vector<uint64_t> single_read;
	string single_fingerprint;
	uint64_t key;
	string fingerprintLeftside;
	string fingerprintRightside;
	vector<int> ATaccumulate;
	int ATsum;
	int len; //length of read

	for (int whichReads = start; whichReads < end; whichReads++) 
	{
		const string& x = Reads[whichReads];
		len = x.length();
		howManyAnAlignmentUnitsInARead = (len - lengthOfAlignmentUnit) / DBA + 1;
		ATaccumulate.resize(len);
		ATsum = 0;
		for (int i = 0; i < len; i++)
		{
			if (x[i] == 'A' || x[i] == 'a' || x[i] == 'T' || x[i] == 't')
				ATsum = ATsum + 1;
			ATaccumulate[i] = ATsum;
		}

		int alignmentUnit_order = 0;
		int translocated_position = 0;
		single_read.resize(howManyAnAlignmentUnitsInARead);


		//______________________________The first type of alignment unit, stored in the file separated by commas________________________________
		while ((translocated_position + lengthOfAlignmentUnit) <= len)
		{
			single_fingerprint = get01Version(ATaccumulate, translocated_position, howManyWindowInAnAlignmentUnit);

			key = generatingKey(single_fingerprint);
			single_read[alignmentUnit_order] = key;

			translocated_position = translocated_position + DBA;
			alignmentUnit_order = alignmentUnit_order + 1;
		}

		v_e[whichReads] = single_read;  // Process by read unit, each read length is 2 * step
		v_e_ID[whichReads] = ReadsID[whichReads];

		//______________________________The second type of alignment unit, stored in the file separated by commas_________________________________

		alignmentUnit_order = 0;
		single_read.resize(2 * DBA_2);

		alignmentUnit_order = 0;
		while (alignmentUnit_order <= (DBA_2 - 1)) // consecutive
		{
			// Calculate the left DBA_2 information fingerprints
			fingerprintLeftside = get01Version(ref(ATaccumulate), alignmentUnit_order, howManyWindowInAnAlignmentUnit);
			key = generatingKey(fingerprintLeftside);

			single_read[2 * alignmentUnit_order] = key;


			// Calculate the right DBA_2 information fingerprints
			fingerprintRightside = get01Version(ref(ATaccumulate), (-1) * alignmentUnit_order - lengthOfAlignmentUnit, howManyWindowInAnAlignmentUnit);
			key = generatingKey(fingerprintRightside);

			single_read[2 * alignmentUnit_order + 1] = key;
			alignmentUnit_order = alignmentUnit_order + 1;

		}

		v_c[whichReads] = single_read; // Process by read unit, each read length is 2 * step
		v_c_ID[whichReads] = ReadsID[whichReads]; 
	}

}


void MAIN_step3(string dir_name, int groupSize_step3_step4, int lastGroupSize_step3_step4, 
	int groups_step3_step4, string grouped_reads_folder_name,
	string grouped_aligningUnit1_folder_name, string grouped_aligningUnit2_folder_name)
{
	mutex mtx;  // Mutex lock
	unsigned int recommanded_threads = std::thread::hardware_concurrency();
	if (recommanded_threads == 0) 
	{
		std::cout << "Unable to retrieve hardware parallelism information" << std::endl;
	}
	else 
	{
		std::cout << "Number of available threads:" << recommanded_threads << std::endl;
	}
	const int thread_s = int(recommanded_threads); // Convert a variable type to a run-time constant (not yet a compile-time constant)


	if (std::filesystem::create_directory(dir_name + "/" + grouped_aligningUnit1_folder_name))
	{
		std::cout << "Successfully created " << dir_name + "/" + grouped_aligningUnit1_folder_name << std::endl;
	}
	else
	{
		std::cerr << "Creation failed (possibly already exists)" << std::endl;
	}


	if (std::filesystem::create_directory(dir_name + "/" + grouped_aligningUnit2_folder_name))
	{
		std::cout << "Successfully created " << dir_name + "/" + grouped_aligningUnit2_folder_name << std::endl;
	}
	else
	{
		std::cerr << "Creation failed (possibly already exists)" << std::endl;
	}


	// "FASTA_lines" is used to record the total number of lines of FASTA file.
	string line;
	vector<string> Reads;
	vector<string> ReadsID;
	vector<vector<uint64_t>> v1_e(groupSize_step3_step4);
	vector<string>v1_e_ID(groupSize_step3_step4);
	vector<vector<uint64_t>> v1_c(groupSize_step3_step4);
	vector<string> v1_c_ID(groupSize_step3_step4);
	

	int loc1; // outer order

	//______________________________________Start grouping and begin the loop_____________________________________

	ifstream inputFile2;

	int i = 0; // current group number
	int u = 0;
	int SIZE; // current group size
	string grouped_reads_file_name;
	string grouped_aligningUnit1_file_name;
	string grouped_aligningUnit2_file_name;

	SIZE = groupSize_step3_step4;
	Reads.resize(SIZE);
	ReadsID.resize(SIZE);

	ofstream outfile1(grouped_reads_file_name);
	ofstream outfile2;
	ofstream outfile3;

	// Store the start and end indices for each thread (left-closed right-open interval)
	std::vector<int> starts, ends;
	starts.reserve(thread_s);
	ends.reserve(thread_s);

	int base = SIZE / thread_s;
	int remainder = SIZE % thread_s;

	int current_start = 0;
	int current_end;
	for (int i = 0; i < thread_s; ++i) 
	{
		// Calculate the end position for the current thread
		current_end = current_start + base + (i < remainder ? 1 : 0);

		starts.push_back(current_start);
		ends.push_back(current_end);

		// Update the starting position for the next thread
		current_start = current_end;
	}

	for (int i = 0; i < groups_step3_step4 * 2 ; i++)
	{
		grouped_reads_file_name = dir_name + "/" + grouped_reads_folder_name + "/" + to_string(i) + ".txt";
		inputFile2.open(grouped_reads_file_name);
		u = 0;
		while (getline(inputFile2, line))
		{
			if (u % 2 == 1)
				Reads[u / 2] = line;
			else
				ReadsID[u / 2] = line;
			u++;
		}
		inputFile2.close();
		//cout << "groups:" << i << endl;

		grouped_aligningUnit1_file_name = dir_name + "/" + grouped_aligningUnit1_folder_name + "/" + to_string(i) + ".txt";
		grouped_aligningUnit2_file_name = dir_name + "/" + grouped_aligningUnit2_folder_name + "/" + to_string(i) + ".txt";

		outfile2.open(grouped_aligningUnit1_file_name);
		outfile3.open(grouped_aligningUnit2_file_name);

		if (i == groups_step3_step4 - 1 || i == groups_step3_step4)
		{
			v1_e.resize(lastGroupSize_step3_step4);
			v1_e_ID.resize(lastGroupSize_step3_step4);
			v1_c.resize(lastGroupSize_step3_step4);
			v1_c_ID.resize(lastGroupSize_step3_step4);

			vector<thread> t;
			int a = SIZE / thread_s;
			for (int i = 0; i < thread_s; i++)
			{
				int start = starts[i];
				int end = ends[i];
				t.emplace_back(gettingAligningUnitAndItsLocation,
					ref(Reads),
					ref(ReadsID),
					start,
					end,
					ref(v1_e),
					ref(v1_e_ID),
					ref(v1_c),
					ref(v1_c_ID)
				);
			}
			for (auto& ti : t)
			{
				ti.join();
			}
		}

		else
		{
			//The sizes of the "Reads" and "ReadsID" are both groupSize
			v1_e.resize(groupSize_step3_step4);
			v1_e_ID.resize(groupSize_step3_step4);
			v1_c.resize(groupSize_step3_step4);
			v1_c_ID.resize(groupSize_step3_step4);

			vector<thread> t;
			int a = SIZE / thread_s;
			for (int i = 0; i < thread_s; i++)
			{
				int start = starts[i];
				int end = ends[i];
				t.emplace_back(gettingAligningUnitAndItsLocation,
					ref(Reads),
					ref(ReadsID),
					start,
					end,
					ref(v1_e),
					ref(v1_e_ID),
					ref(v1_c),
					ref(v1_c_ID)
				);
			}
			for (auto& ti : t)
			{
				ti.join();
			}
		}

		// write in outfile2, outfile3
		loc1 = 0;
		for (vector<uint64_t> x : v1_e)
		{
			for (uint64_t xx : x)
			{
				outfile2 << xx << ",";
			}
			outfile2 << endl;
			outfile2 << v1_e_ID[loc1] << endl;
			loc1 = loc1 + 1;
		}

		loc1 = 0;
		for (vector<uint64_t> x : v1_c)
		{
			for (uint64_t xx : x)
			{
				outfile3 << xx << ",";
			}
			outfile3 << endl;
			outfile3 << v1_c_ID[loc1] << endl;
			loc1 = loc1 + 1;
		}

		outfile2.close();
		outfile3.close();
		
		if (i + 1 == groups_step3_step4 - 1 || i + 1 == groups_step3_step4) // if next group is lastGroup
		{
			SIZE = lastGroupSize_step3_step4;
			Reads.resize(SIZE);
			ReadsID.resize(SIZE);

			base = SIZE / thread_s;
			remainder = SIZE % thread_s;

				starts.clear();
				ends.clear();
				starts.reserve(thread_s);
				ends.reserve(thread_s);
				current_start = 0;
				for (int i = 0; i < thread_s; ++i)
				{
					current_end = current_start + base + (i < remainder ? 1 : 0);

					starts.push_back(current_start);
					ends.push_back(current_end);

					current_start = current_end;
				}
			}
		else
		{
				SIZE = groupSize_step3_step4;
				Reads.resize(SIZE);
				ReadsID.resize(SIZE);

				base = SIZE / thread_s;
				remainder = SIZE % thread_s;

				starts.clear();
				ends.clear();
				starts.reserve(thread_s);
				ends.reserve(thread_s);
				current_start = 0;
				for (int i = 0; i < thread_s; ++i)
				{
					current_end = current_start + base + (i < remainder ? 1 : 0);

					starts.push_back(current_start);
					ends.push_back(current_end);

					current_start = current_end;
				}
		}
	}
}