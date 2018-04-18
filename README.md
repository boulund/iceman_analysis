# Analysis of Iceman samples
Analysis of samples from the Iceman.

Snakemake workflow designed to be used with a reference background sample, in
this case a lung tissue sample (1018). The reference background sample is 
(unfortunately) hardcoded in some of the workflow steps, because of some Snakemake
quirk I don't understand. The primary sample, a gut sample, is sample 224.

# Dependencies
The workflow uses several third-party tools and databases, and expects all of
these to be available in the running environment (i.e. to be in `$PATH`).

## Tools
* BBMap (and several other BBTools):
 - `bbduk.sh`, quality filtering, adapter trimming, etc.
 - `bbcountunique.sh`, assess kmer saturation
 - `bbmerge.sh`, read merging
 - `bbmap.sh`, read mapping
 - `tadpole.sh`, assembly
 - `translate6frames.sh`, six-frame translation
 - `pileup.sh`, coverage statistics
* HMMER3, for TIGRFAM annotation
* PMDTools, for filtering ancient reads
* MetaPhlAn2, for taxonomic profiling
* MicrobeCensus, for average genome size estimation
* mapDamage, assess ancient DNA indications

Some additional custom tools are available in the `scripts` directory, all of
which are written for Python 3.6+, and some require standard scientific packages 
such as Pandas and Matplotlib. 

## Databases
* MEGARes, database of antibiotic resistance genes
* TIGRFAMs, database of protein families 

# Download and run
Get the workflow by cloning this repository:

```
git clone --recursive https://github.com/boulund/iceman_analysis.git
```

This will produce a folder called `iceman_analysis` in your current directory. Edit `config.yaml` so that
all paths are correct for your system/environment. Then run the pipeline:

```
snakemake --snakefile /path/to/iceman_analysis/Snakefile --configfile /path/to/iceman_analysis/config.yaml --jobs N
```

Set the number of jobs (`N`) to the number of available CPU cores to maximize parallelization potential. 



# Tunapuco_Matses
In the folder `tunapuco_matses` is a self-contained workflow for running
similar analyses of another dataset intended for comparison with the Iceman
data. 

To run the Tunapuco_Matses workflow, configure relevant paths and the input
filename pattern in the `config.yaml` for that workflow, then run:

```
snakemake --use-conda --jobs N
```

to run the entire workflow for all samples found in the defined input folder. The workflow
uses [conda](www.conda.io) to make the workflow easier to reproduce.
