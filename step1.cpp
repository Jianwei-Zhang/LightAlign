#include"f.h"
#include"globals.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <string>

using namespace std;
namespace fs = filesystem;

string getFileExtension(const std::string& fileName) 
{
    size_t dotPos = fileName.find_last_of('.');
    if (dotPos != std::string::npos && dotPos != 0) 
    {
        return fileName.substr(dotPos + 1);
    }
    return "";
}

void MAIN_step1(string dir_name, string inputfile_step1, string outputfile_step1) 
{
    const int groupSize = 10000;
    vector<string> block;
    string line;
    string sequenceName;

    vector<pair<string, int>>sorting;
    unordered_map<string, uint64_t> name_bias;

    //===========fasta or fastq=============
    string extension = getFileExtension(inputfile_step1);

    // mkdir
    try
    {
        if (fs::create_directory(dir_name))
            cout << "Directory created successfully!" << endl;
        else
            cout << "Directory already exists or failed to create." << endl;
    }
    catch (const exception& e)
    {
        cerr << "Error creating directory: " << e.what() << endl;
        return;
    }
    ifstream inputFile(inputfile_step1);
    string sortedReads = dir_name + "/" + "sortedFASTA.fasta";
    ofstream outputFile(sortedReads);
    // Count the total number of lines in a FASTA/FASTQ file 
    int total = 0;

    if (extension == "fastq"|| extension == "fq")
    {
        long long sum_len = 0;

        total = 0;
        string tempLine;
        string nameline;
        while (getline(inputFile, tempLine))
        {
            if ( total % 4 == 0)
            {
                name_bias[tempLine] = static_cast<uint64_t>(inputFile.tellg());
                nameline = tempLine;
            }
            else if( total % 4 == 1)
                sorting.push_back({ nameline , size(tempLine) });

            total++;
        }
        inputFile.clear();
        inputFile.seekg(0);

        // Sort in descending order by uint64_t. This step consumes significant memory when the number of sequences is large.
        sort(sorting.begin(), sorting.end(),
            [](const auto& a, const auto& b) 
            {
                return a.second > b.second;
            });

        cout <<  total / 4 << " reads." << endl;

        // Output the original sequencing file sorted by length
        for (auto i : sorting)
        {
            outputFile << i.first << endl;
            inputFile.seekg(name_bias[i.first]);
            getline(inputFile, line);
            outputFile << line << endl;
        }
        inputFile.close();
        outputFile.close();

        total = total / 2;
    }
    
    else if (extension == "fasta"|| extension == "fa")
    {
        long long sum_len = 0;

        total = 0;
        string tempLine;
        string nameline;
        while (getline(inputFile, tempLine))
        {
            if (total % 2 == 0)
            {
                name_bias[tempLine] = static_cast<uint64_t>(inputFile.tellg()); // Line offset
                nameline = tempLine;
            }
            else
                sorting.push_back({ nameline , size(tempLine) });
            total++;
        }
        inputFile.clear();
        inputFile.seekg(0);

        sort(sorting.begin(), sorting.end(),
            [](const auto& a, const auto& b)
            {
                return a.second > b.second;
            });

        cout << total / 2 << " reads." << endl;

        for (auto i : sorting)
        {
            outputFile << i.first << endl;
            inputFile.seekg(name_bias[i.first]);
            getline(inputFile, line);
            outputFile << line << endl;
        }
        inputFile.close();
        outputFile.close();
    }
    
    else
    {
        cout << "only fasta and fastq." << endl;
    }
    
     inputFile.open(sortedReads);
     outputFile.open(outputfile_step1);


     const int groups = total / groupSize;
     const int lastGroupSize = total % groupSize;
     const int least_overlap = lengthOfAlignmentUnit + 2 * DBA;

     for (int i = 0; i <= groups; i++)
     {
         const int currentGroupSize = (i == groups) ? lastGroupSize : groupSize;
         block.resize(currentGroupSize);

         for (int lineNum = 0; lineNum < currentGroupSize; lineNum++)
         {
             if (!getline(inputFile, line))
             {
                 break;
             }
             block[lineNum] = line;
         }

         for (int loc = 0; loc < currentGroupSize; loc++)
         {
             const string& x = block[loc];
             const int len = x.length();

             if (loc % 2 == 1 && len >= least_overlap)
             {
                 sequenceName = ">" + block[loc - 1].substr(1);
                 outputFile << sequenceName << endl;
                 outputFile << x << endl;
             }
         }

         if (i != groups) 
         {
             for (int skip = 0; skip < groupSize - currentGroupSize; skip++)
             {
                 if (!getline(inputFile, line)) break;
             }
         }
     }

     inputFile.close();
     outputFile.close();

     fs::path sorted = dir_name + "/" + "sortedFASTA.fasta";
     fs::remove(sorted);
 }