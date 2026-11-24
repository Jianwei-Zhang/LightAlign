LightAlign Linux source and build instructions

Build from the repository root:
    bash linux/build_linux.sh

Run from the repository root with 12 worker threads:
    ./linux/bin/LightAlign -i reads.fasta -O output_directory \
        -g 10200 -w 30 -e 0.02 -d 87 -l 1200 -t 12

Requirements: a Linux C++17 compiler (g++), pthreads and C++ filesystem support.
The -l 1200 example is for the revised eukaryotic experiments; the prokaryotic
reruns used -l 840. Use a fresh output directory for every run because the
overlap output uses append mode. FASTA/FASTQ input must be uncompressed.

This linux/src directory contains the server Linux source. Windows source and
build files elsewhere in the repository describe a separate build.
The binary archive contains the existing Linux x86_64 executable, its SHA256
checksum and runtime dependency information. A rebuild may have a different hash.
