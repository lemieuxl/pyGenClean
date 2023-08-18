# Marker call rate

The _marker_call_rate_ QC module consists of a single [script](#main-script). It
uses _Plink_ to remove markers with low call rate.

Use the following command to access the multiple scripts of the
_marker_call_rate_ QC module.

```shell-session
$ pyGenClean marker-call-rate --help
usage: pyGenClean marker-call-rate [-h] [-v] {run} ...

Remove markers with poor call rate.

options:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

subcommands:
  Below is a list of tools in the marker-call-rate module. Note that 'run'
  executes the main marker-call-rate pipeline.

  {run}
    run          Remove markers with poor call rate.
```

## Main script

The main script (accessed using the `run` subcommand) uses _Plink_
to remove heterozygous haploid genotypes.

```shell-session
$ pyGenClean marker-call-rate run --help
usage: pyGenClean marker-call-rate run [-h] --bfile FILE [--geno FLOAT]
                                       [--plink-1.07] [--out FILE]

Remove markers with poor call rate.

options:
  -h, --help    show this help message and exit

Input File:
  --bfile FILE  The input file prefix (will find the plink binary files by
                appending the prefix to the .bim, .bed and .fam files,
                respectively).

Options:
  --geno FLOAT  The missingness threshold (remove SNPs with more than x
                percent missing genotypes). [Default: 0.02]
  --plink-1.07  Use original Plink (version 1.07)

Output File:
  --out FILE    The prefix of the output files. [default: clean_geno]
```

### Input Files

This module uses _Plink_'s binary file format (`bed`, `bim` and `fam` files) for
the source data set (the data of interest).

### Procedure

Here are the steps performed by the module:

1. Run Plink with the `geno` option.
2. Compare the two `bim` files (before and after the Plink `geno`
   analysis).

### Output files

Here is a comprehensive list of all the possible output files for each of the
steps described above.

!!! note
    The output prefix shown is the one by default (_i.e._ `clean_geno`).

`clean_geno.{bed,bim,fam}`
: The dataset with markers having a high missing rate removed (according to a
  user defined threshold).

`clean_geno.removed_snps`
: The list of markers that have a high missing rate (above a user defined
  threshold).
