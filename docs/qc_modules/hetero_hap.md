# Heterozygous haploid

The _hetero_hap_ QC module consists of a single [script](#main-script). It uses
_Plink_ to remove heterozygous haploid genotypes.

Use the following command to access the multiple scripts of the _hetero_hap_ QC
module.

```shell-session
$ pyGenClean hetero-hap --help
usage: pyGenClean hetero-hap [-h] [-v] {run} ...

Remove heterozygous haploid (males on chromosome X).

options:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

subcommands:
  Below is a list of tools in the hetero-hap module. Note that 'run'
  executes the main hetero-hap pipeline.

  {run}
    run          Remove heterozygous haploid (males on chromosome X).
```

## Main script

The main script (accessed using the `run` subcommand) uses _Plink_
to remove heterozygous haploid genotypes.

```shell-session
$ pyGenClean hetero-hap run --help
usage: pyGenClean hetero-hap run [-h] --bfile FILE [--plink-1.07] [--out FILE]

Remove heterozygous haploid (males on chromosome X).

options:
  -h, --help    show this help message and exit

Input File:
  --bfile FILE  The input file prefix (will find the plink binary files by
                appending the prefix to the .bim, .bed and .fam files,
                respectively).

Options:
  --plink-1.07  Use original Plink (version 1.07)

Output File:
  --out FILE    The prefix of the output files. [default:
                without_hh_genotypes]
```

### Input Files

This module uses _Plink_'s binary file format (`bed`, `bim` and `fam` files) for
the source data set (the data of interest).

### Procedure

Here are the steps performed by the module:

1. Remove the heterozygous haploid genotypes uing _Plink_.

### Output files

Here is a comprehensive list of all the possible output files for each of the
steps described above.

!!! note
    The output prefix shown is the one by default (_i.e._
    `without_hh_genotypes`).

`without_hh_genotypes.{bed,bim,fam}`
: The data set with heterozygous haploid genotypes removed.
