# Energy Efficient Quantum Simulations
## Contents

This repository contains the experiments that were ran to collect the data for the CIUK poster – *Energy Efficient Quantum Simulations*. It has [QuEST](https://github.com/QuEST-Kit/QuEST) and [ITensor](https://github.com/ITensor/ITensor) repositories linked, as those are necessary to run the circuit. The implemented circuits are the *Quantum Fourier Transform (QFT)* and a *Random circuit (RAND)* from the [Google's quantum advantage paper](https://www.nature.com/articles/s41586-019-1666-5). 

The implementations are stored in `quest-projects` and `itensor-projects` directories. ITensor also has some additional helper and I/O functions stored in the `itensor-projects/helpers` directory. The `jobs` directory contains scripts for running the experiments on ARCHER2. The results are scraped from the output files and processed in the `results` directory. 


## Building and running

First, set up QuEST and ITensor according to the instructions in their original repositories. ITensor needs to be pre-compiled, while QuEST is directly linked when compiling the input circuit. To compile QuEST implementations, run 

`cd quest-projects && ./build.sh circuit-name`, 

where `circuit-name` can be either `qft` or `rand`. To compile ITensor implementations, run 

`cd itensor-projects/{circuit-name} && make`, 

where `circuit-name` can be either `qft` or `bench` (which stands for the random circuit). 

All executives will be built in the local `bin` directory (i.e. `itensor-projects/bin` or `quest-projects/bin`). 

To run the programs via SLURM, use the appropriate script from the `jobs` directory. Make sure that the appropriate output directory has been created in the same location as the script. 
