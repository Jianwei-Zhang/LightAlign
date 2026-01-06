#include"f.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <cstdio>

// Check if sequences in a FASTA file are in multiline format (only checking the first few lines)
bool has_multiline_sequences(const std::string& input_file, int max_sequences_to_check) {
    std::ifstream input(input_file, ios::binary);
    if (!input) {
        std::cerr << "Error: Could not open file: " << input_file << std::endl;
        return false;
    }

    std::string line;
    bool in_sequence = false;
    int sequence_count = 0;
    int lines_in_current_sequence = 0;

    while (std::getline(input, line) && sequence_count < max_sequences_to_check) {
        // skip the empty line
        if (line.empty()) continue;

        if (line[0] == '>') {
            if (in_sequence && lines_in_current_sequence > 1) {
                input.close();
                return true;
            }

            sequence_count++;
            in_sequence = true;
            lines_in_current_sequence = 0;
        }
        else if (in_sequence) {
            lines_in_current_sequence++;
            if (lines_in_current_sequence > 1) {
                input.close();
                return true;
            }
        }
    }

    input.close();
    return false;
}

// Format FASTA file (merge multiline sequences into single lines)
void format_fasta_to_single_line(const std::string& input_file) {
    std::string temp_file = input_file + ".tmp";

    std::ifstream input(input_file, ios::binary);
    if (!input) {
        std::cerr << "Error: Cannot open input file: " << input_file << std::endl;
        return;
    }

    std::ofstream output(temp_file);
    if (!output) {
        std::cerr << "Error: Cannot create temporary file" << std::endl;
        input.close();
        return;
    }

    std::string line, sequence;
    bool in_sequence = false;

    while (std::getline(input, line)) {
        if (!line.empty() && line[0] == '>') {
            if (in_sequence) {
                output << sequence << '\n';
                sequence.clear();
            }
            output << line << '\n';
            in_sequence = true;
        }
        else if (in_sequence) {
            for (char c : line) {
                if (!std::isspace(static_cast<unsigned char>(c))) {
                    sequence += c;
                }
            }
        }
    }

    if (in_sequence && !sequence.empty()) {
        output << sequence << '\n';
    }

    input.close();
    output.close();

    // replace original file with the tmp file
    std::remove(input_file.c_str());
    std::rename(temp_file.c_str(), input_file.c_str());
}

void prepro(const std::string& input_file, int check_sequences) {
    if (has_multiline_sequences(input_file, check_sequences)) {
        std::cout << "Format FASTA file: " << input_file << std::endl;
        format_fasta_to_single_line(input_file);
    }
    else {
        std::cout << "The FASTA file is already in single-line format: " << input_file << std::endl;
    }
}