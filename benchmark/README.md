## DeepREx: Deep learning-based predictor of Residue EXposure

### Benchmarking of DeepREx

This section contains source and instructions on how to reproduce the DeepREx
results on a blind test set comprising 200 protein sequences and presented
in the following paper:

Manfredi et al. (2021) DeepREx: Deep learning-based prediction of Residue
solvent Exposure from protein sequence (submitted)

To speed up the process, results will be generated using pre-computed
multiple sequence alignments and HMM files using the HHBlits program and the
the [UniClust30 database](http://wwwuser.gwdg.de/~compbiol/uniclust/2020_02/),
version 2020_02.

To regenerate HHblits files for each protein, download the
[UniClust30 database](http://wwwuser.gwdg.de/~compbiol/uniclust/2020_02/)
and run:

```
$ hhblits -i a_protein.fa -d UniRef30_2020_02 -n 2 -oa3m a_protein.a3m -ohhm a_protein.hhm -o /dev/null
```

#### Step 1: clone the DeepREx repository

Clone the repository using Git running:

```
$ git clone https://github.com/BolognaBiocomp/deeprex.git
```

#### Step 2: run benchmarking code

After cloning, move to the benchmark folder and run the benchmarking script.
Scoring measures will be printed to standard output after completion.

```
$ cd deeprex/benchmark
$ ./run_benchmark.sh
SENSITIVITY:    0.828137666740364
PRECISION:      0.8016979987871438
F1:     0.8147033767740298
Q2:     0.8181155036828809
MCC:    0.6365228220061144
```

The script will also produce a file named bts_deeprex.tsv containing all
predictions and ground truth values for the 200 proteins in the benchmark.
