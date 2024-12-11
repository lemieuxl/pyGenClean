# Related samples

The _related_samples_ QC module consists of two scrips: a
[main](#main-script) to find related samples, and single on one to include and
randomly selected samples from related clusters.

Use the following command to access the multiple scripts of the
_related_samples_ QC module.

```shell-session
$ pyGenClean related-samples --help
usage: pyGenClean related-samples [-h] [-v] {run,merge-related-samples} ...

Finds related samples according to IBS values.

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

subcommands:
  Below is a list of tools in the related-samples module. Note that 'run'
  executes the main related-samples pipeline.

  {run,merge-related-samples}
    run                 Finds related samples according to IBS values.
    merge-related-samples
                        Merges related samples according to IBS.
```

## Main script

The main script (accessed using the `run` subcommand) uses _Plink_
to find samples possilby related.

```shell-session
$ pyGenClean related-samples run --help
usage: pyGenClean related-samples run [-h] --bfile FILE [--genome-only]
                                      [--min-nb-snp INT]
                                      [--indep-pairwise STR STR STR]
                                      [--maf FLOAT] [--ibs2-ratio FLOAT]
                                      [--nb-threads N] [--plink-1.07]
                                      [--out FILE]

Finds related samples according to IBS values.

options:
  -h, --help            show this help message and exit

Input files:
  --bfile FILE          The input file prefix (will find the plink binary
                        files by appending the prefix to the .bim, .bed and
                        .fam files, respectively.)

Options:
  --genome-only         Only create the genome file
  --min-nb-snp INT      The minimum number of markers needed to compute IBS
                        values. [10000]
  --indep-pairwise STR STR STR
                        Three numbers: window size, window shift and the r2
                        threshold. [['50', '5', '0.1']]
  --maf FLOAT           Restrict to SNPs with MAF >= threshold. [0.05]
  --ibs2-ratio FLOAT    The initial IBS2* ratio (the minimum value to show in
                        the plot. [0.8]
  --nb-threads N        The number of threads for this analysis (no effect
                        when using plink 1.07). [1]
  --plink-1.07          Use original Plink (version 1.07). Note that this will
                        be slow, as Plink 1.07 doesn't support multi
                        threading.

Output files:
  --out FILE            The prefix of the output files. [ibs]
```

### Input Files

This module uses _Plink_'s binary file format (`bed`, `bim` and `fam` files) for
the source data set (the data of interest). It also requires the plate
organization for each samples.
