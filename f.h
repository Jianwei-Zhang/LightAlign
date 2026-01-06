#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <regex>
#include <string>
#include <cctype>
#include <climits>
#include<mutex>
#include <unordered_set>
#include <atomic>
#include <queue>
#include <future>
#include <sstream>
#include"globals.h"

using namespace std;

bool has_multiline_sequences(const std::string& input_file, int max_sequences_to_check=3);
void format_fasta_to_single_line(const std::string& input_file);
void prepro(const std::string& input_file, int check_sequences=3);

void printHelp();

//____________________________________________fastq_into_fastq____________________________________________
void MAIN_step1(string dir_name, string inputfile_step1, string outputfile_step1);

//______________________________________________allReadsAndItsID__________________________________________
string reverseComplement(string DNAsequence);

tuple<int, int, int>  MAIN_step2(string dir_name, string inputfile_step2, int groupSize_step3_step4);

uint64_t generatingKey(const std::string& AligningUnit);

//__________________________________________aligningUnitAndItsLocation______________________________________
string get01Version(vector<int>& ATaccumulate, int startpoint, int howManyWindowInAnAlignmentUnit);// Calculate the alignmentUnit corresponding to a single reads fragment

void gettingAligningUnitAndItsLocation(
	vector<string>& Reads,
	vector<string>& ReadsID,
	int start,
	int end,
	vector<vector<uint64_t>>& v_e, vector<string>& v_e_ID,
	vector<vector<uint64_t>>& v_c, vector<string>& v_c_ID);

void MAIN_step3(string dir_name, int groupSize_step3_step4, int lastGroupSize_step3_step4,
	int groups_step3_step4, string grouped_reads_folder_name,
	string grouped_aligningUnit1_folder_name, string grouped_aligningUnit2_folder_name);

//______________________________________________NeedlemanWunsch_________________________________________
tuple<string, string, int> needlemanWunsch
(
	const string& seqA,
	const string& seqB,
	int8_t match,
	int8_t mismatch,
	int8_t gap
);

//______________________________________________overlap___________________________________________________

bool contains(const std::string& mainStr, const std::string& subStr);

tuple<int, int, int, string> difference(const string& str1, const string& str2);

void calculateBytesOfBlocks(const std::string& file, // inputPath
	vector<long long>& groupSizes_bytes, // output
	int groupSize,
	int groups,
	int lastGroupSize);

std::string extractID(const string& str);

void seperate(
	vector<uint64_t>& row,
	vector <string>& blockxy_AligningUnit1AndItsLocation,
	vector<vector<uint64_t>>& keyxy_evenAligningUnit,
	vector<string>& blockxy_location_of_evenAligningUnit,
	int groupSize);

void seperate2(
	vector<uint64_t>& row,
	vector <string>& blockxy_AligningUnit2AndItsLocation,
	vector<vector<uint64_t>>& keyxy_consecutiveAligningUnit,
	vector<string>& blockxy_location_of_consecutiveAligningUnit,
	int groupSize);

void calculate_partitions(
	int group_size,
	int num_threads,
	std::vector<int>& start_inside1,
	std::vector<int>& end_inside1
);

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
	string outfile_filted,
	bool is_lastGroup
);

void alignment(string filtedFile,
	ofstream& outfile_overlap,
	ofstream& outfile_score,
	unordered_map<string, long long>& reads_hashtable,
	string outputfile_step1
);


string trim(const string& str);


long long countLinesFast(const std::string& filename);

void read_in_data_for_allReads(ifstream& inputFile,
	string& dir_name, string& grouped_reads_folder_name,
	vector<std::string>& block11_allReads,
	vector<std::string>& block12_allReads,
	vector<int>& which_group_had_alignment_inside,
	int groupSize, int firstGroup,
	int groupSize_step3_step4, int groups_step3_step4, int lastGroupSize_step3_step4);

void read_in_data_for_all_AligningUnitAndItsLocation(ifstream& inputFile,
	string& dir_name, string& grouped_aligningUnit_folder_name,
	vector<std::string>& block11_AligningUnitAndItsLocation,
	vector<std::string>& block12_AligningUnitAndItsLocation,
	vector<int>& which_group_had_alignment_inside,
	int groupSize, int firstGroup,
	int groupSize_step3_step4, int groups_step3_step4, int lastGroupSize_step3_step4);

void read_in_data_for_allReads_secondGroup(ifstream& inputFile,
	string& dir_name, string& grouped_reads_folder_name,
	vector<string>& block21_allReads,
	vector<string>& block22_allReads,
	int groupSize, int secondGroup,
	int groupSize_step3_step4, int groups_step3_step4, int lastGroupSize_step3_step4);

void read_in_data_for_all_AligningUnitAndItsLocation_secondGroup(ifstream& inputFile,
	string& dir_name, string& grouped_aligningUnit_folder_name,
	vector<string>& block21_AligningUnitAndItsLocation,
	vector<string>& block22_AligningUnitAndItsLocation,
	int groupSize, int secondGroup,
	int groupSize_step3_step4, int groups_step3_step4, int lastGroupSize_step3_step4);