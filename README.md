# HamSampler

HamSampler is a C program that samples RNA sequences from a given RNA secondary structure with a given Boltzmann distribution as well as with Hamming distance filtration. It takes a secondary structure S, a reference sequence t_0 and a given distance d as input, and generates RNA sequences t that are Boltzmann distributed having a fixed Hamming distance to t_0. 

The program is written in C. Please cite the software as specified at the bottom of the paper. 


### Prerequisites

GNU Automake

C Standard Library


### Installing

Use the command: 

```
./configure
make
```
The executable file will be "Sampler" in ./src. 


## Running 

### Running HamSampler 

There are several options for HamSampler. Use the following command to see all options. 
```
./src/Sampler -h  
```

To specify an input file, use the command
```
./src/Sampler -i input.in 
```
The default input file is "input.in". 

To specify an output file, use the command
```
./src/Sampler -o output.out 
```
The default output file is "output.out". 

To set the number of sampled sequences to be c (a positive integer), use the command
```
./src/Sampler -m c
```
The default number of sampled sequences is c = 1000. 

To set the Hamming distance to the reference sequence to be h (a positive integer smaller than the sequence length), use the command
```
./src/Sampler -d h 
```
The default distance is h = 5. 

To compute the EMR-spectrum up to distance h, use the command
```
./src/Sampler -i input.in -o output.out -r -d h 
```
This option is disabled by default. The user needs to specify the option "-r" to enable this functionality. Then the EMR value at distance h of the input structure can be found in "output.out". The program also provides the EMR value of randomly sampled sequences that are compatible with the given structure. 

### Input file style

An example of an input can be found in "input.in". The input file consists of two lines. The first line is a secondary structure in dot/bracket form. We recommend the length of the input structure to be below 500. Here is an example of a valid input secondary structure.  
```
.(((((((.((((((......))))))((((((.....)))))).....((((((.....)))))))))))))....
```

The second line is the reference sequence. It is not necessary for it to be compatible with the input structure. 
Here is an example of a valid reference sequence for the above input structure.  
```
UGGCCCCCAACCCGCCGUACCGCGGGUUGCCGACACCGUCGGCGUGGUACCCCGCCCUACGCGGGGGGGGGCCCCAC
```

Currently we do not support input with multiple structures. 

### Output file

If no EMR-spectrum (no "-r" option) is set, the user can find c entries of sampled sequences in the output file. We also provide their frequencies after each sampled sequence. The included output file "output.out" is an example by running the commands
```
./src/Sampler -d 10 
```

If the EMR-spectrum (with "-r" option) is set, the sampled sequences are not displayed. The user can find the EMR value at distance h in the output file. The included output file "spectrum.out" is an example by running the commands
```
./src/Sampler -o spectrum.out -r -d 10 
```


## Contact

If you have any question or bug reports about HamSampler, please feel free to conatct the author Fenix Huang at fenixprotoss@gmail.com.  




