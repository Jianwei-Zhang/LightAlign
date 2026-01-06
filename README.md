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
```LightAlign.exe -O [Output path for results and intermediate files] -i [Input file path, FASTA/FASTQ]```

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
## Notes:
LightAlign executes in 5 steps. Upon completion, "All steps done." will be displayed. The final output is the PAF file in the output path.
In tests conducted so far, LightAlign consumes the same amount of memory and has roughly the same runtime on both Linux and Windows.  

## Limitations
Currently, LightAlign is limited to the assembly of HiFi reads.
Despite its extremely low RAW usage, LightAlign is slower than current tools like minimap2.
## Contact
For any questions or suggestions, please reach out to liujian_HZAU@outlook.com.
