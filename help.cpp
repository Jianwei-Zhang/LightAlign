#include"f.h"
#include <iostream>
#include <string>

using namespace std;

void printHelp() {
    cout << "Usage: LightAlign.exe -O [output path] -i [input path]\n\n"
        << "Options:\n"
        << "-O      requires a directory path.\n"
        << "-i      requires an input file path.\n"
        << "-w requires an integer value, the recommended value is 25~35.\n"
        << "-l requires an integer value, the recommended value is 800 ~ 900 for prokaryotes and 1000 ~ 1100 for eukaryotes.\n"
        << "-e requires a float value, the recommended value is 0.02.\n"
        << "-g requires a int value, the recommended value is 10200.\n"
        << "-d requires a int value, the recommended value is 70~90, depends on the length of reads.\n"
        << " -p  preprocess fasta file\n";
}