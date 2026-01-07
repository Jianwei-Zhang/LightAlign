# LightAlign
LightAlign is a memory-efficient alignment tool for HiFi data that utilizes sequence fuzzy feature and reduces the peak memory usage during overlaps detection. It outputs Pairwise Alignment Format (PAF), when combined with miniasm, LightAlign enables high-quality bacterial genome assembly with memory usage below 1 GB, while also making eukaryotic genome assembly feasible on standard personal computers or even laptops.  
LightAlign is fully compliant with the C++17 language standard.  
## Installation
### windows system:
To use LightAlign, download all the source code to your computer, place it in the same project, compile this project using a C++ IDE such as Visual Studio, and generate the executable (.exe) file to start using the tool; or you can also download LightAlign.exe to your computer.
### Linux system:
One Compilation Example Using GCC 12.2.0:  
```
g++ -std=c++17 \  
-isystem /public/home/software/opt/bio/software/GCC/12.2.0/include/c++/12.2.0 \  
-ILightAlign/include \  
LightAlign/src/help.cpp \  
LightAlign/src/step4_5.cpp \  
LightAlign/src/step3.cpp \  
LightAlign/src/step2.cpp \  
LightAlign/src/step1.cpp \  
LightAlign/src/functions.cpp \  
LightAlign/src/MurmurHash3.cpp \  
LightAlign/src/preprocess_fasta.cpp \  
-o LightAlign_executable \  
-pthread \  
-lstdc++fs  
```
## Tutorial
### LightAlign Basic Usage:  
```LightAlign.exe -i [Input file path, FASTA/FASTQ] -O [Output path for results and intermediate files]```

### LightAlign Detailed Usage

LightAlign accepts the following command-line options:  

| Option | Type    | Description | Default |
|--------|---------|-------------|---------|
| `-h`   | flag    | Show help documentation and exit. | |
| `-i`   | string  | Path to the input file (required). | |
| `-O`   | string  | Path for the output file (required). | |
| `-w`   | INT     | Window size for alignment. | 30 |
| `-l`   | INT     | Alignment unit length / minimum overlap length. Recommended: 800–900 for prokaryotes, 1000–1100 for eukaryotes. | 810 |
| `-e`   | FLOAT   | Maximum allowed error rate for alignments. | 0.02 |
| `-g`   | INT     | Number of reads per group in the group alignment phase. Reducing this value lowers memory usage (current memory usage is consistently <1 GB for prokaryotic datasets). Keep default unless ultra-low memory is required. | 10200 |
| `-d`   | INT     | DBA (Dynamic Bandwidth Adjustment). Adjust to around 80 when average HiFi read length is <10 kb. | 87 |  
## Quick Test Example

You can validate LightAlign using the public dataset **SRR32655838**(PacBio HiFi data of Thalassospira sp. MIT1370) from NCBI SRA ([link](https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR32655838)). The following command demonstrates a typical alignment operation and its expected outcome within a pipeline.

### Run Command & Expected Results
1.  **Run LightAlign** on the dataset with the following command. It typically uses about **0.54 GB** of memory.
    ```bash
    LightAlign.exe -i SRR32655838/SRR32655838.fasta -O [Output path] -d 75 -l 900
    ```
2.  **Use the output for assembly** with tools like **Miniasm**. The following command shows how the results can be piped into the next step:
    ```bash
    miniasm -f SRR32655838.fasta SRR32655838.paf > output.gfa
    ```
3.  **Final Assembly Metric**: When processing this dataset with the above LightAlign parameters followed by Miniasm, the contig **N50 in the resulting GFA file is typically around 3.11 Mb**.  (Parameters `-d 75` and `-l 900` are optimized for this specific dataset. Actual memory usage and N50 may vary slightly depending on the system and tool versions.)

## Notes:
LightAlign executes in 5 steps. Upon completion, "All steps done." will be displayed. The final output is the PAF file in the output path.
In tests conducted so far, LightAlign consumes the same amount of memory and has roughly the same runtime on both Linux and Windows.  

## Limitations
Currently, LightAlign is limited to the assembly of HiFi reads.
## Contact
For any questions or suggestions, please reach out to liujian_HZAU@outlook.com.
