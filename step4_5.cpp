#include"f.h"
#include<iostream>
#include<fstream>
#include<string>
#include <thread> // Multithreading; Gets the number of threads the program can use
#include<sstream> // separate strings by commas
#include<functional>
#include<vector>
#include<algorithm>//the find function
#include <memory>
#include <mutex>
#include <time.h>
#include<filesystem>
#include<cmath>// The square root function
#include<cstdlib> // Accepts parameters entered on the command line
#include"globals.h"
#include <unordered_set>
#include<unordered_map> // creating hash table
#include <omp.h>
#include <immintrin.h>
#include<set>
#include <cctype>

using namespace std;

//default value
int DBA = 87;
int DBA_2 = DBA;
int step = 30;
int lengthOfAlignmentUnit = 810; // Prokaryotic and eukaryotic datasets require different lengths
int window = 30;
double error_rate = 0.02;
bool cigar = true; // Whether to output the CIGAR string; speed is basically the same

void createFile(const string& filePath)
{
	std::ofstream file(filePath);
}

int main(int argc, char* argv[])
{
	// Check for user input
	string dir_name;
	string inputfile_step1;
	int groupSize_step3_step4 = 10200;

	// Parse command line arguments
	//bool preprocess_flag = false; // By default, multi-line FASTA files do not require preprocessing. Use the -p parameter if preprocessing is needed.
	bool show_help = false;
	for (int i = 1; i < argc; ++i) {
		std::string arg = argv[i];

		if (argc == 1) {
			printHelp();
			return 0;
		}

		if (arg == "-h" || arg == "--help") {
			show_help = true;
			break;
		}

		else if (arg == "-O") {
			if (i + 1 < argc) {
				dir_name = argv[++i];
			}
			else {
				std::cerr << "-O requires a directory path." << std::endl;
				return 1;
			}
		}

		else if (arg == "-i") {
			if (i + 1 < argc) {
				inputfile_step1 = argv[++i];
			}
			else {
				std::cerr << "-i requires an input file path." << std::endl;
				return 1;
			}
		}

		else if (arg == "-w") {
			if (i + 1 < argc) {
				window = std::atoi(argv[++i]);
			}
			else {
				std::cerr << "-w requires an integer value, the recommended value is 25~35." << std::endl;
				return 1;
			}
		}

		else if (arg == "-l") {
			if (i + 1 < argc) {
				lengthOfAlignmentUnit = std::atoi(argv[++i]);
			}
			else {
				std::cerr << "-l requires an integer value, the recommended value is 800 ~ 900 for prokaryotes and 1000 ~ 1100 for eukaryotes." << std::endl;
				return 1;
			}
		}

		else if (arg == "-e") {
			if (i + 1 < argc) {
				error_rate = std::atof(argv[++i]);
			}
			else {
				std::cerr << "-e requires a float value, the recommended value is 0.02." << std::endl;
				return 1;
			}
		}

		else if (arg == "-g") {
			if (i + 1 < argc) {
				groupSize_step3_step4 = std::atof(argv[++i]);
			}
			else {
				std::cerr << "-g requires a int value, the recommended value is 10200." << std::endl;
				return 1;
			}
		}

		else if (arg == "-d") {
			if (i + 1 < argc) {
				DBA = std::atof(argv[++i]);
			}
			else {
				std::cerr << "-d requires a int value, the recommended value is 70~90, depends on the length of reads." << std::endl;
				return 1;
			}
		}
		/*
		else if (arg == "-p") {
			preprocess_flag = true;
		}
		*/
	}
	if (show_help) {
		printHelp();
		return 0;
	}
	/*
	if (preprocess_flag) 
	{
		prepro(inputfile_step1);
	}
	*/
	if (dir_name.empty() || inputfile_step1.empty()) 
	{
		std::cerr << "Error: -O and -i parameters are required." << std::endl;
		std::cerr << "Usage: " << argv[0] << " -O <directory> -i <input_file> [-w <value>] [-l <value>] [-e <value>] [-d]" << std::endl;
		return 1;
	}

	std::cout << "DBA: " << DBA << std::endl;
	std::cout << "step: " << step << std::endl;
	std::cout << "lengthOfAlignmentUnit: " << lengthOfAlignmentUnit << std::endl;
	std::cout << "window: " << window << std::endl;
	std::cout << "error_rate: " << error_rate << std::endl;

	// step1
	prepro(inputfile_step1);

	string outputfile_step1 = dir_name + "/" + "filteredFASTA.fasta";
	cout << "============Step1:Input Standardization===========" << endl;
	MAIN_step1(dir_name, inputfile_step1, outputfile_step1);

	// step2
	string inputfile_step2 = outputfile_step1;
	cout << "==========Step2:Complementary Strand Generation==========" << endl;
	auto [FASTA_lines, lastGroupSize_step3_step4, groups_step3_step4] = MAIN_step2(dir_name, inputfile_step2, groupSize_step3_step4);

	int howManyReads = (groups_step3_step4 - 1) * groupSize_step3_step4 + lastGroupSize_step3_step4; // Not including complementary chains, the actual number of reads should be twice of it
	// step3
	string grouped_reads_folder_name = "grouped_reads";
	string grouped_aligningUnit1_folder_name = "grouped_aligningUnit1";
	string grouped_aligningUnit2_folder_name = "grouped_aligningUnit2";

	cout << "============Step3:Fuzzy Feature Extraction============" << endl;
	MAIN_step3(dir_name, groupSize_step3_step4, lastGroupSize_step3_step4, groups_step3_step4, grouped_reads_folder_name, grouped_aligningUnit1_folder_name, grouped_aligningUnit2_folder_name);

	cout << "============Step4:Correlated Read Pair Screening============" << endl;
	mutex fileMutex;
	vector<int> which_group_had_alignment_inside; //firstGroup and secondGroup that have overlap_inside

	int windowNumberOfAlignmentUnit = (lengthOfAlignmentUnit - window) / step + 1;
	int temp;

	unsigned int recommanded_threads = std::thread::hardware_concurrency();
	if (recommanded_threads == 0)
	{
		std::cout << "Unable to retrieve hardware parallelism information" << std::endl;
	}
	else
	{
		std::cout << "Number of available threads: " << recommanded_threads << std::endl;
	}
	int thread_s = int(recommanded_threads); // Convert a variable type to a run-time constant (not yet a compile-time constant)

	for (int i = 0; i < thread_s; i++)
	{
		string filtedFile = dir_name + "/filted_" + to_string(i) + ".txt";
		createFile(filtedFile);
	}


	string line;
	//_____________________________________________Open file_______________________________________________

	ifstream inputFile4;
	std::cout << "groups: " << groups_step3_step4 << endl;
	std::cout << "lastGroupSize: " << lastGroupSize_step3_step4 << endl;


	vector <string>block11_AligningUnit1AndItsLocation(2 * groupSize_step3_step4); // One read contains two lines
	vector <string>block12_AligningUnit1AndItsLocation(2 * groupSize_step3_step4);
	vector <string>block21_AligningUnit1AndItsLocation(2 * groupSize_step3_step4);
	vector <string>block22_AligningUnit1AndItsLocation(2 * groupSize_step3_step4);
	ifstream inputFile2;


	vector <string>block11_AligningUnit2AndItsLocation(2 * groupSize_step3_step4);
	vector <string>block12_AligningUnit2AndItsLocation(2 * groupSize_step3_step4);
	vector <string>block21_AligningUnit2AndItsLocation(2 * groupSize_step3_step4);
	vector <string>block22_AligningUnit2AndItsLocation(2 * groupSize_step3_step4);
	ifstream inputFile3;

	//_____________________________________________prepare before loop_______________________________________________

	int i = 0;
	int loc = 0; // which line

	// Split the string into evenAligningUnit
	vector<vector<uint64_t>> key11_evenAligningUnit(groupSize_step3_step4);
	vector<string> block11_location_of_evenAligningUnit(groupSize_step3_step4);
	vector<vector<uint64_t>> key12_evenAligningUnit(groupSize_step3_step4);
	vector<string> block12_location_of_evenAligningUnit(groupSize_step3_step4);// first dimension
	vector<vector<uint64_t>> key21_evenAligningUnit(groupSize_step3_step4);
	vector<string> block21_location_of_evenAligningUnit(groupSize_step3_step4);// first dimension
	vector<vector<uint64_t>> key22_evenAligningUnit(groupSize_step3_step4);
	vector<string> block22_location_of_evenAligningUnit(groupSize_step3_step4);// first dimension

	vector<uint64_t> row_in_a1_former;
	vector<vector<string>> row2D_in_a1_loc_former;
	vector<uint64_t> row_in_a1_latter;
	vector<vector<string>> row2D_in_a1_loc_latter;
	vector<uint64_t> row_in_a2_former(2 * DBA_2);
	vector<vector<string>> row2D_in_a2_loc_former(DBA_2 * 2, vector<string>(2));
	vector<uint64_t> row_in_a2_latter(2 * DBA_2);
	vector<vector<string>> row2D_in_a2_loc_latter(DBA_2 * 2, vector<string>(2));


	// Split the string into consecutiveAligningUnit
	vector<vector<uint64_t>> key11_consecutiveAligningUnit(groupSize_step3_step4);
	vector<vector<uint64_t>> key12_consecutiveAligningUnit(groupSize_step3_step4);
	vector<vector<uint64_t>> key21_consecutiveAligningUnit(groupSize_step3_step4);
	vector<string> block21_location_of_consecutiveAligningUnit(groupSize_step3_step4);
	vector<vector<uint64_t>> key22_consecutiveAligningUnit(groupSize_step3_step4);
	vector<string> block22_location_of_consecutiveAligningUnit(groupSize_step3_step4);

	string processed_l;

	vector<uint64_t> roww(2 * DBA_2); // store a line of struct

	vector<vector<vector<string>>> location_of_evenAligningUnit;
	vector<vector<vector<string>>> location_of_consecutiveAligningUnit;

	// The inside1 parameter, not required for every loop
	vector<int> start_inside1(thread_s);
	vector<int> end_inside1(thread_s);
	// The inside2 parameter, not required for every loop
	vector<int> start_inside2(thread_s);
	vector<int> end_inside2(thread_s);
	// between parameter
	vector<int> start_between(thread_s);
	vector<int> end_between(thread_s);

	int length_of_side; // Base length Initial value
	int h; // The sum of the first n h's
	float area; // Each small trapezoid area

	bool is_lastGroup = false; // Whether it has cycled to the last group

	unordered_set<string> used_Reads;

	//___________________________________________what lastGroup need_____________________________________________

	vector <string> all_AligningUnit1AndItsLocation;
	vector <string> all_AligningUnit2AndItsLocation;
	// Split the string into alignment units
	istringstream one_Read;
	int a_or_l = 0; // When read as a alignment unit, it determines that it is aligningUnit or location
	istringstream is;
	string buffer;
	string l;
	clock_t start2, end2;

	//________________________________________start clocking______________________________________________
	start2 = clock();

	//___________________________________________Start the first loop_____________________________________________
	for (int firstGroup = 0; firstGroup < groups_step3_step4; firstGroup++)
	{
		cout << "First group is " << firstGroup << " Now" << endl;

		// Multiple types of data are simultaneously read and processed as in-memory data types
		read_in_data_for_all_AligningUnitAndItsLocation(ref(inputFile2),
			ref(dir_name), ref(grouped_aligningUnit1_folder_name),
			ref(block11_AligningUnit1AndItsLocation), ref(block12_AligningUnit1AndItsLocation),
			ref(which_group_had_alignment_inside),
			groupSize_step3_step4, firstGroup,
			groupSize_step3_step4, groups_step3_step4, lastGroupSize_step3_step4);

		inputFile2.close();
		inputFile3.close();
		calculate_partitions(groupSize_step3_step4, thread_s, ref(start_inside1), ref(end_inside1));
		calculate_partitions(groupSize_step3_step4, thread_s, ref(start_inside2), ref(end_inside2));


		// read blockxy_AligningUnit1AndItsLocation into blockxy_evenAligningUnit, blockxy_location_of_evenAligningUnit
		seperate(ref(row_in_a1_former),
			ref(block11_AligningUnit1AndItsLocation),
			ref(key11_evenAligningUnit),
			ref(block11_location_of_evenAligningUnit),
			groupSize_step3_step4);


		//================================Generating hash table===================================
		unordered_map<uint64_t, vector<vector<string> > > hash_table;
		int hashSize = DBA * 3 * groupSize_step3_step4;
		int aligningUnit_order_firstRead;
		hash_table.reserve(hashSize);

		for (int x = 0; x < groupSize_step3_step4; x++) // x represents the index of firstRead in the current block11_allReads
		{
			aligningUnit_order_firstRead = 0;
			for (uint64_t key : key11_evenAligningUnit[x])
			{
				auto it = hash_table.find(key);  // checking if the key already exists
				if (it != hash_table.end())
				{
					hash_table[key].push_back({ to_string(x),
						to_string(aligningUnit_order_firstRead * DBA) });
					// Each vector<string> contains: {index of firstRead in the current block11_allReads, start position of the alignment unit}
				}
				else
				{
					// Stores the starting position of the alignment unit on current first read
					hash_table.emplace(key, vector<vector<string>>{{to_string(x),
						to_string(aligningUnit_order_firstRead* DBA)}});
				}
				aligningUnit_order_firstRead = aligningUnit_order_firstRead + 1;
			}
		}


		//________________________________________________second loops________________________________________________
		for (int secondGroup = firstGroup; secondGroup < groups_step3_step4; secondGroup++) // firstGroup + 1
		{
			// Whether it has cycled to the last group
			if (secondGroup == groups_step3_step4 - 1)
			{
				calculate_partitions(lastGroupSize_step3_step4, thread_s, ref(start_inside2), ref(end_inside2));

				// Resize the data block based on the size of the last group
				block21_AligningUnit2AndItsLocation.resize(2 * lastGroupSize_step3_step4);
				block22_AligningUnit2AndItsLocation.resize(2 * lastGroupSize_step3_step4);
				key21_consecutiveAligningUnit.resize(lastGroupSize_step3_step4);
				block21_location_of_consecutiveAligningUnit.resize(lastGroupSize_step3_step4);
				key22_consecutiveAligningUnit.resize(lastGroupSize_step3_step4);
				block22_location_of_consecutiveAligningUnit.resize(lastGroupSize_step3_step4);

				//______________________________________read in data for second group______________________________________

				// Multiple types of data are simultaneously read and processed as in-memory data types

				read_in_data_for_all_AligningUnitAndItsLocation_secondGroup(ref(inputFile3),
					ref(dir_name), ref(grouped_aligningUnit2_folder_name),
					ref(block21_AligningUnit2AndItsLocation), ref(block22_AligningUnit2AndItsLocation),
					lastGroupSize_step3_step4, secondGroup,
					groupSize_step3_step4, groups_step3_step4, lastGroupSize_step3_step4);

				// read in AligningUnit1AndItsLocation and AligningUnit2AndItsLocation
				thread d[2];
				d[0] = thread(seperate2,
					ref(row_in_a2_former),
					ref(block21_AligningUnit2AndItsLocation),
					ref(key21_consecutiveAligningUnit),
					ref(block21_location_of_consecutiveAligningUnit),
					lastGroupSize_step3_step4);
				d[1] = thread(seperate2,
					ref(row_in_a2_latter),
					ref(block22_AligningUnit2AndItsLocation),
					ref(key22_consecutiveAligningUnit),
					ref(block22_location_of_consecutiveAligningUnit),
					lastGroupSize_step3_step4);

				for (int i = 0; i < 2; i++)
				{
					d[i].join();
				}
			}

			else
			{
				// start_inside2 is related to secondGroups
				calculate_partitions(groupSize_step3_step4, thread_s, ref(start_inside2), ref(end_inside2));
				//______________________________________read in data for second group______________________________________
				// Multiple types of data are simultaneously read and processed as in-memory data types
				read_in_data_for_all_AligningUnitAndItsLocation_secondGroup(
					ref(inputFile3),
					ref(dir_name),
					ref(grouped_aligningUnit2_folder_name),
					ref(block21_AligningUnit2AndItsLocation),
					ref(block22_AligningUnit2AndItsLocation),
					groupSize_step3_step4, secondGroup,
					groupSize_step3_step4, groups_step3_step4, lastGroupSize_step3_step4);


				inputFile2.close(); // Close the file after reading
				inputFile3.close(); // Close the file after reading

				// read in AligningUnit1AndItsLocation and AligningUnit2AndItsLocation
				thread d[2];
				d[0] = thread(seperate2,
					ref(row_in_a2_former),
					ref(block21_AligningUnit2AndItsLocation),
					ref(key21_consecutiveAligningUnit),
					ref(block21_location_of_consecutiveAligningUnit),
					groupSize_step3_step4);
				d[1] = thread(seperate2,
					ref(row_in_a2_latter),
					ref(block22_AligningUnit2AndItsLocation),
					ref(key22_consecutiveAligningUnit),
					ref(block22_location_of_consecutiveAligningUnit),
					groupSize_step3_step4);

				for (int i = 0; i < 2; i++)
				{
					d[i].join();
				}
			}


			//  Whether it has cycled to the last group
			if (secondGroup == groups_step3_step4 - 1)
			{
				is_lastGroup = true;
				for (int i = 0; i < thread_s; i++) // range of second_Read
				{
					start_between[i] = (lastGroupSize_step3_step4 / thread_s) * i;
					end_between[i] = (lastGroupSize_step3_step4 / thread_s) * (i + 1);
				}
			}
			else
			{
				is_lastGroup = false;
				for (int i = 0; i < thread_s; i++)
				{
					start_between[i] = (groupSize_step3_step4 / thread_s) * i;
					end_between[i] = (groupSize_step3_step4 / thread_s) * (i + 1);
				}
			}

			// overlap_between
			vector<std::thread> t;
			for (int i = 0; i < thread_s; i++) //thread_s
			{
				string filtedFile = dir_name + "/filted_" + to_string(i) + ".txt";
				t.emplace_back(overlap_between,
					ref(block11_location_of_evenAligningUnit),
					ref(block21_location_of_consecutiveAligningUnit),
					ref(block22_location_of_consecutiveAligningUnit),
					ref(key21_consecutiveAligningUnit),
					ref(key22_consecutiveAligningUnit),
					ref(hash_table),
					howManyReads,
					firstGroup, secondGroup, groups_step3_step4,
					start_between[i],
					end_between[i],
					windowNumberOfAlignmentUnit,
					groupSize_step3_step4,
					lastGroupSize_step3_step4,
					filtedFile,
					is_lastGroup
				);
			}
			for (auto& ti : t)
			{
				ti.join();
			}


			if (is_lastGroup == true) // Reset the data block size
			{
				block21_AligningUnit1AndItsLocation.resize(2 * groupSize_step3_step4);
				block22_AligningUnit1AndItsLocation.resize(2 * groupSize_step3_step4);
				block21_AligningUnit2AndItsLocation.resize(2 * groupSize_step3_step4);
				block22_AligningUnit2AndItsLocation.resize(2 * groupSize_step3_step4);

				key21_evenAligningUnit.resize(groupSize_step3_step4);
				block21_location_of_evenAligningUnit.resize(groupSize_step3_step4);
				key22_evenAligningUnit.resize(groupSize_step3_step4);
				block22_location_of_evenAligningUnit.resize(groupSize_step3_step4);
				key21_consecutiveAligningUnit.resize(groupSize_step3_step4);
				block21_location_of_consecutiveAligningUnit.resize(groupSize_step3_step4);
				key22_consecutiveAligningUnit.resize(groupSize_step3_step4);
				block22_location_of_consecutiveAligningUnit.resize(groupSize_step3_step4);
			}
			end2 = clock(); // time out
			//cout << "time = " << double(end2 - start2) / CLOCKS_PER_SEC << "s" << endl;
		}

		// Adding complementary strand groups to the hash table for each group achieves the overlap_inside effect

		hash_table.clear();
	}

	vector<string>().swap(block11_AligningUnit1AndItsLocation);
	vector<string>().swap(block12_AligningUnit1AndItsLocation);
	vector<string>().swap(block21_AligningUnit1AndItsLocation);
	vector<string>().swap(block22_AligningUnit1AndItsLocation);

	vector<string>().swap(block11_AligningUnit2AndItsLocation);
	vector<string>().swap(block12_AligningUnit2AndItsLocation);
	vector<string>().swap(block21_AligningUnit2AndItsLocation);
	vector<string>().swap(block22_AligningUnit2AndItsLocation);

	vector<vector<uint64_t>>().swap(key11_evenAligningUnit);
	block11_location_of_evenAligningUnit.clear();// first dimension
	vector<vector<uint64_t>>().swap(key12_evenAligningUnit);
	block12_location_of_evenAligningUnit.clear();// first dimension
	vector<vector<uint64_t>>().swap(key21_evenAligningUnit);
	block21_location_of_evenAligningUnit.clear();// first dimension
	vector<vector<uint64_t>>().swap(key22_evenAligningUnit);
	block22_location_of_evenAligningUnit.clear();// first dimension

	vector<vector<uint64_t>>().swap(key11_consecutiveAligningUnit);
	vector<vector<uint64_t>>().swap(key12_consecutiveAligningUnit);
	vector<vector<uint64_t>>().swap(key21_consecutiveAligningUnit);
	vector<vector<uint64_t>>().swap(key22_consecutiveAligningUnit);
	vector<string>().swap(block21_location_of_consecutiveAligningUnit);
	vector<string>().swap(block22_location_of_consecutiveAligningUnit);

	//_________________________Delete intermediate files_________________________
	filesystem::path dir_path;
	dir_path = dir_name + "/" + grouped_reads_folder_name;
	filesystem::remove_all(dir_path);

	dir_path = dir_name + "/" + grouped_aligningUnit1_folder_name;
	filesystem::remove_all(dir_path);

	dir_path = dir_name + "/" + grouped_aligningUnit2_folder_name;
	filesystem::remove_all(dir_path);

	cout << "============Step5:Segmented Banded Alignment============" << endl;

	unordered_map<string, long long> reads_hashtable;
	reads_hashtable.reserve(howManyReads * 2); // including complementary chains
	vector<long long> index;

	//_____________________________________Read the line offset file______________________________________
	ifstream inputFile1(dir_name + "/" + "index.txt", ios::binary);
	while (getline(inputFile1, line))
	{
		index.push_back(stoll(line));
	}
	inputFile1.close();

	// Here we need to read each sequence name line and its corresponding line offset into a hash table.
	inputFile1.open(inputfile_step2, ios::binary);
	i = 0;
	while (getline(inputFile1, line))
	{
		if (i % 2 == 0) {
			string trimmed_line = trim(line);
			reads_hashtable.emplace(trimmed_line, index[i + 1]);
		}
		i++;
	}
	inputFile1.close();

	string outputfile = dir_name + "/overlap.paf";
	ofstream outfile_overlap(outputfile, ios::app);
	outputfile = dir_name + "/score.txt";
	ofstream outfile_score(outputfile, ios::app);

	vector<thread> t;
	vector<string> temp_overlap_files;
	vector<string> temp_score_files;
	for (int i = 0; i < thread_s; i++)
	{
		string filtedFile = dir_name + "/filted_" + to_string(i) + ".txt";
		t.emplace_back(alignment,
			filtedFile,
			ref(outfile_overlap),
			ref(outfile_score),
			ref(reads_hashtable),
			inputfile_step2
		);
	}
	for (auto& ti : t)
	{
		ti.join();
	}

	// remove filted_i.txt
	for (int i = 0; i < thread_s; i++)
	{
		string filtedFile = dir_name + "/filted_" + to_string(i) + ".txt";
		dir_path = filtedFile;
		filesystem::remove_all(dir_path);
	}

	cout << "All steps done." << endl;
	system("pause");

	return 0;
}
