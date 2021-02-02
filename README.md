## DeepREx: Deep learning-based predictor of Residue EXposure

### The DeepREx Docker image

Image availbale on DockerHub [https://hub.docker.com/r/bolognabiocomp/deeprex](https://hub.docker.com/r/bolognabiocomp/deeprex)

#### Usage of the image

The first step to run DeepREx Docker container is the pull the container image. To do so, run:

```
$ docker pull bolognabiocomp/deeprex
```

Now the DeepREx Docker image is installed in your local Docker environment and ready to be used. To show DeepREx help page run:

```
$ docker run bolognabiocomp/deeprex -h

usage: deeprex.py [-h] -f FASTA -d HHBLITS_DB -o OUTF [-a CPUS]

DeepREx: Deep learning-based predictor of Residue EXposure

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        The input multi-FASTA file name
  -d HHBLITS_DB, --hhblits-database HHBLITS_DB
                        The HHBlits database file
  -o OUTF, --outf OUTF  The output TSV file
```
The program accepts three mandatory arguments:
- The path of the input FASTA file containing protein sequences to be predicted;
- The hhblits database prefix;
- The output file where predictions will be stored.

Let's now try a concrete example. First of all, create a new working folder and cd into it:

```
$ mkdir my_deeprex_example
$ cd my_deeprex_example
```

Then, let's download an example sequence from UniProtKB, e.g. the TUridylate kinase protein 52 from Prochlorococcus marinus with accession :

```
$ wget https://www.uniprot.org/uniprot/Q46GQ4.fasta
```
We also need a valid HHBlits database. In our DeepREx server we use the [UniClust30 database](http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/). Here, for simplicity, we'll use a smaller database e.g. the [Pfam HHblits database](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_33.1.tar.gz).
Let's download and unpack the database:

```
$ wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_33.1.tar.gz
$ tar xvzf pfamA_33.1.tar.gz
```

Now, we are ready to predict the solvent exposure of our input protein. Run:

```
$ docker run -v $(pwd):/data/ bolognabiocomp/deeprex -f Q46GQ4.fasta -d pfam -o Q46GQ4.tsv
```

In the example above, we are mapping the current program working directory ($(pwd)) to the /data/ folder inside the container. This will allow the container to see the external FASTA file Q46GQ4.fasta and the HHblits database.
The file Q46GQ4.tsv now contains the DeepREx prediction, in TSV format:
```
$ cat Q46GQ4.tsv
sp|Q46GQ4|PYRH_PROMT	1	M	Exposed	0.74
sp|Q46GQ4|PYRH_PROMT	2	T	Exposed	0.63
sp|Q46GQ4|PYRH_PROMT	3	Y	Exposed	0.59
sp|Q46GQ4|PYRH_PROMT	4	S	Buried	0.02
sp|Q46GQ4|PYRH_PROMT	5	R	Buried	0.06
sp|Q46GQ4|PYRH_PROMT	6	A	Buried	0.96
sp|Q46GQ4|PYRH_PROMT	7	L	Buried	0.98
sp|Q46GQ4|PYRH_PROMT	8	I	Buried	0.99
sp|Q46GQ4|PYRH_PROMT	9	K	Buried	0.92
sp|Q46GQ4|PYRH_PROMT	10	L	Buried	0.99
sp|Q46GQ4|PYRH_PROMT	11	S	Buried	0.66
sp|Q46GQ4|PYRH_PROMT	12	G	Buried	0.73
sp|Q46GQ4|PYRH_PROMT	13	E	Exposed	0.39
sp|Q46GQ4|PYRH_PROMT	14	A	Buried	0.61
sp|Q46GQ4|PYRH_PROMT	15	L	Buried	0.9
sp|Q46GQ4|PYRH_PROMT	16	M	Exposed	0.22

```
Columns are as follows:
- Column 1: the protein ID/accession as reported in the FASTA input file;
- Column 2: the residue position;
- Column 3: the residue name;
- Column 4: predicted solvent exposure: Buried or Exposed;
- Column 5: reliability index associated to the prediction.

### Install and use DeepREx from source

Source code available on GitHub at [https://github.com/BolognaBiocomp/deeprex](https://github.com/BolognaBiocomp/deeprex)

#### Installation and configuration

DeepREx is designed to run on Unix/Linux platforms. The software was written using the Python programming language and it was tested under the Python version 3.

To obtain DeepREx, clone the repository from GitHub:

```
$ git clone https://github.com/BolognaBiocomp/deeprex.git
```

This will produce a directory “deeprex”. Before running deeprex you need to set and export a variable named DEEPREX_ROOT to point to the deeprex installation dir:
```
$ export DEEPREX_ROOT='/path/to/deeprex'
```

Before running the program, you need to install DeepREx dependencies. We suggest to use Conda (we suggest [Miniconda3](https://docs.conda.io/en/latest/miniconda.html)) create a Python virtual environment and activate it.

To create a conda env for deeprex:

```
$ conda create -n deeprex
```
To activate the environment:

```
$ conda activate deeprex
```

The following dependencies are required:

- biopython (version 1.78 tested)
- Keras (version 2.4.3 tested)
- HHSuite (version 3.3.0)

To install all requirements:

```
$ conda install nomkl keras biopython
$ conda install -c conda-forge -c bioconda hhsuite
```

Now you are able to use deeprex(see next Section). Remember to keep the environment active.
If you whish, you can copy the “deeprex.py” script to a directory in the users' PATH.

#### Usage

The program accepts three mandatory arguments:

- The path of the input FASTA file containing protein sequences to be predicted;
- The hhblits database prefix;
- The output file where predictions will be stored.

We also need a valid HHBlits database. In our DeepREx server we use the [UniClust30 database](http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/). Here, for simplicity, we'll use a smaller database e.g. the [Pfam HHblits database](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_33.1.tar.gz).
Let's download and unpack the database:

```
$ wget http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/pfamA_33.1.tar.gz
$ tar xvzf pfamA_33.1.tar.gz
```

As an example, run the program on the example FASTA file contained in the folder "example":

```
$ ./deepsig.py -f example/Q46GQ4.fasta -d pfam -o example/Q46GQ4.tsv
```
Once the prediction is done, the TSV output should look like the following:

```
$ cat example/Q46GQ4.tsv
sp|Q46GQ4|PYRH_PROMT	1	M	Exposed	0.74
sp|Q46GQ4|PYRH_PROMT	2	T	Exposed	0.63
sp|Q46GQ4|PYRH_PROMT	3	Y	Exposed	0.59
sp|Q46GQ4|PYRH_PROMT	4	S	Buried	0.02
sp|Q46GQ4|PYRH_PROMT	5	R	Buried	0.06
sp|Q46GQ4|PYRH_PROMT	6	A	Buried	0.96
sp|Q46GQ4|PYRH_PROMT	7	L	Buried	0.98
sp|Q46GQ4|PYRH_PROMT	8	I	Buried	0.99
sp|Q46GQ4|PYRH_PROMT	9	K	Buried	0.92
sp|Q46GQ4|PYRH_PROMT	10	L	Buried	0.99
sp|Q46GQ4|PYRH_PROMT	11	S	Buried	0.66
sp|Q46GQ4|PYRH_PROMT	12	G	Buried	0.73
sp|Q46GQ4|PYRH_PROMT	13	E	Exposed	0.39
sp|Q46GQ4|PYRH_PROMT	14	A	Buried	0.61
sp|Q46GQ4|PYRH_PROMT	15	L	Buried	0.9
sp|Q46GQ4|PYRH_PROMT	16	M	Exposed	0.22
....
```
Columns are as follows:
- Column 1: the protein ID/accession as reported in the FASTA input file;
- Column 2: the residue position;
- Column 3: the residue name;
- Column 4: predicted solvent exposure: Buried or Exposed;
- Column 5: reliability index associated to the prediction.

Please, reports bugs to: castrense.savojardo2@unibo.it
