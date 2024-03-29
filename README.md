# RICseqlib

## Quick start
Clone this repository with the following command:

```{bash}
git clone https://github.com/malyshev-andrey/RICseqlib.git .
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

