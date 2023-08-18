# Flag MAF

The _flag_maf_ QC module consists of a single [script](#main-script). It uses
_Plink_ to flag monomorphic markers (_i.e._ a minor allele frequency (MAF) of
0). Note that no marker gets excluded from this step.

Use the following command to access the multiple scripts of the _flag_maf_ QC
module.

```shell-session
$ pyGenClean flag-maf --help
usage: pyGenClean flag-maf [-h] [-v] {run} ...

Flag markers with MAF of 0 (monomorphic).

options:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

subcommands:
  Below is a list of tools in the flag-maf module. Note that 'run' executes
  the main flag-maf pipeline.

  {run}
    run          Flag markers with MAF of 0 (monomorphic).
```

## Main script

The main script (accessed using the `run` subcommand) uses _Plink_
to find monomorphic markers (_i.e._ markers with a MAF of 0).

```shell-session
$ pyGenClean flag-maf run --help
usage: pyGenClean flag-maf run [-h] --bfile FILE [--plink-1.07] [--out FILE]

Flag markers with MAF of 0 (monomorphic).

options:
  -h, --help    show this help message and exit

Input File:
  --bfile FILE  The input file prefix (will find the plink binary files by
                appending the prefix to the .bim, .bed and .fam files,
                respectively.

Options:
  --plink-1.07  Use original Plink (version 1.07)

Output File:
  --out FILE    The prefix of the output files. [default: flag_maf_0]
```

### Input Files

This module uses _Plink_'s binary file format (`bed`, `bim` and `fam` files) for
the source data set (the data of interest).

### Procedure

Here are the steps performed by the module:

1. Computes the frequencies using Plink.
2. Finds markers with MAF of 0, and saves them in a file.

### Output files

Here is a comprehensive list of all the possible output files for each of the
steps described above.

!!! note
    The output prefix shown is the one by default (_i.e._ `flag_maf_0`).

`flag_maf_0.frq`
: The frequency of each marker in the source dataset.

`flag_maf_0.list`
: The list of markers with a minor allele frequency of zero.
