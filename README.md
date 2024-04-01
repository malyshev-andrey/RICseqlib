# RICseqlib

## Quick start
Clone this repository with the following command:

```{bash}
git clone https://github.com/malyshev-andrey/RICseqlib.git ./RICseqlib
```

## Dependencies

- [Python3](https://github.com/python/cpython)
- [Miniconda](https://docs.anaconda.com/free/miniconda/index.html)
- [Snakemake](https://anaconda.org/bioconda/snakemake)

## Usage

Script main.sbatch can be used to run pipeline via [Slurm](https://github.com/SchedMD/slurm) job scheduler (Snakemake installation is done automatically in this case so
only conda is needed):

```{bash}
sbatch main.sbatch
```

**Alternatively**, direct snakemake start is available (the recommended number of
cores is 12):

```
snakemake --profile profile --cores 12
```

## Configuration

The file `pipeline_config.yml` can be used to specify samples names, links or paths to `.fastq` files with raw RIC-seq reads, a genome assembly version and annotation. Also, it's possible to modify the arguments of programs included in the pipeline (see respective sections of yaml file).
