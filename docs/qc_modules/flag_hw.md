# Flag Hardy-Weinberg

The _flag_hw_ QC module consists of a single [script](#main-script). It uses
_Plink_ to flag markers failing Hardy-Weinberg equilibrium. Note that no marker
gets excluded from this step.

Use the following command to access the multiple scripts of the _flag_hw_ QC
module.

```shell-session
$ pyGenClean flag-hw --help
usage: pyGenClean flag-hw [-h] [-v] {run} ...

Flags markers failing Hardy Weinberg equilibrium.

options:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

subcommands:
  Below is a list of tools in the flag-hw module. Note that 'run' executes
  the main flag-hw pipeline.

  {run}
    run          Flags markers failing Hardy Weinberg equilibrium.
```

## Main script

The main script (accessed using the `run` subcommand) uses _Plink_
to find markers who fail Hardy-Weinberg equilibrium test.

```shell-session
$ pyGenClean flag-hw run --help
usage: pyGenClean flag-hw run [-h] --bfile FILE [--hwe FLOAT] [--plink-1.07]
                              [--out FILE]

Flags markers failing Hardy Weinberg equilibrium.

options:
  -h, --help    show this help message and exit

Input File:
  --bfile FILE  The input file prefix (will find the plink binary files by
                appending the prefix to the .bim, .bed and .fam files,
                respectively).

Options:
  --hwe FLOAT   The Hardy-Weinberg equilibrium threshold. [default: 1e-4]
  --plink-1.07  Use original Plink (version 1.07)

Output File:
  --out FILE    The prefix of the output files. [default: flag_hw]
```

### Input Files

This module uses _Plink_'s binary file format (`bed`, `bim` and `fam` files) for
the source data set (the data of interest).

### Procedure

Here are the steps performed by the module:

1. Compute the number of markers in the input file.
2. Compute the Bonferroni threshold.
3. Run _Plink_ to find failed markers for HW equilibrium with the Bonferroni
   threshold.
4. Run _Plink_ to find failed markers for HW equilibrium with the default
   threshold.
5. Compare the two marker lists (Bonferroni and default threshold) and finds
   markers that are between the two thresholds.

### Output files

The output files of each of the steps described above are as follow (note that
the output prefix shown is the one by default [i.e. flag_hw]):

`flag_hw.threshold_X.Xe-X.{bed,bim,fam}`
: The data set containing only the markers passing the HW equilibrium test
  (above the Bonferroni threshold).

`flag_hw.threshold_1e-4.{bed,bim,fam}`
: The data set containing only the markers passing the HW equilibrium test
  (above the genome wide significance threshold of $1 \times 10^{−4}$). This
  value can be modified at the command line.

`flag_hw.snp_flag_threshold_X.Xe-X`
: The list of markers failing the HW equilibrium test for the Bonferroni
  threshold.

`flag_hw.snp_flag_threshold_1e-4`
: The list of markers failing the HW equilibrium test for the genome wide
  significance threshold of $1 \times 10^{−4}$. This value can be modified at
  the command line.

`flag_hw.snp_flag_threshold_between_1e-4-X.Xe-X`
: The list of markers failing the HW equilibrium test at a threshold between the
  Bonferroni and the genome wide significance thresholds, so that you can
  exclude only the ones that have a lower p value than the Bonferroni threshold.
