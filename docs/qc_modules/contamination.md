# Contamination

The _contamination_ QC module consists of a single [main script](#main-script).
It tries to find contaminated samples using the
[BAFRegress](https://genome.sph.umich.edu/wiki/BAFRegress) tool.

Use the following command to access the script of the _contamination_ QC module.

```shell-session
$ pyGenClean contamination --help
usage: pyGenClean contamination [-h] [-v] {run} ...

Check for sample contamination using BAFRegress.

options:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

subcommands:
  Below is a list of tools in the contamination module. Note that 'run'
  executes the main contamination pipeline.

  {run}
    run          Check for sample contamination using BAFRegress.
```

## Main script

The main script (accessed using the `run` subcommand) uses
[BAFRegress](https://genome.sph.umich.edu/wiki/BAFRegress) to estimate sample's
contamination.

```shell-session
$ pyGenClean contamination run --help
usage: pyGenClean contamination run [-h] --bfile FILE --raw-dir DIR
                                    [--colsample COL] [--colmarker COL]
                                    [--colbaf COL] [--colab1 COL]
                                    [--colab2 COL]
                                    [--estimate-threshold FLOAT] [--nb-cpu NB]
                                    [--max-samples-per-job NB] [--plink-1.07]
                                    [--out FILE]

Check for sample contamination using BAFRegress.

options:
  -h, --help            show this help message and exit

Input File:
  --bfile FILE          The input file prefix (will find the plink binary
                        files by appending the prefix to the .bim, .bed and
                        .fam files, respectively).

Raw Data:
  --raw-dir DIR         Directory containing the raw data (one file per
                        sample, where the name of the file (minus the
                        extension) is the sample identification number.
  --colsample COL       The sample column. [default: Sample Name]
  --colmarker COL       The marker column. [default: SNP Name]
  --colbaf COL          The B allele frequency column. [default: B Allele
                        Freq]
  --colab1 COL          The AB Allele 1 column. [default: Allele1 - AB]
  --colab2 COL          The AB Allele 2 column. [default: Allele2 - AB]

Options:
  --estimate-threshold FLOAT
                        The estimate threshold for which a sample is
                        considered contaminated. [>0.01]
  --nb-cpu NB           The number of CPU to use. [default: 1]
  --max-samples-per-job NB
                        The maximum number of samples per task. [100]
  --plink-1.07          Use original Plink (version 1.07). Note that this will
                        be slow, as Plink 1.07 doesn't support multi
                        threading.

Output File:
  --out FILE            The prefix of the output files. [default:
                        contamination]
```

### Input files

This module uses Plink's binary file format (`bed`, `bim` and `fam` files) for
the source data set (the data of interest). It also uses intensities file (one
per sample to test, option `--raw-dir`) usually provided by the genotyping
platform.

#### Intensities

[BAFRegress](https://genome.sph.umich.edu/wiki/BAFRegress) uses one intensity
file (_CSV_ format) per sample. This file is generallty generated directly by
the genotyping platform.

The keyword `[Data]` should be present on the line before the _header_ of the
file. Then, there should be one column for each of those required columns (see
table below for a list of required columns, their default values and the option
used to change the default values).

| Required column                 | Default value   | Option        |
| :------------------------------ | :-------------- | :------------ |
| Sample (ID) information         | `Sample Name`   | `--colsample` |
| Marker (name) information       | `SNP Name`      | `--colmarker` |
| B allele frequency column (BAF) | `B Allele Freq` | `--colbaf`    |
| A allele                        | `Allele1 - AB`  | `--colab1`    |
| B allele                        | `Allele2 - AB`  | `--colab2`    |

Below are some lines for a _csv_ report containing only the required columns for
[BAFRegress](https://genome.sph.umich.edu/wiki/BAFRegress).

```text
[Data]
SNP Name,Sample Name,Allele1 - AB,Allele2 - AB,B Allele Freq
1:103380393,Sample_1,B,B,1.0
1:109439680,Sample_1,A,A,0.014306434
1:110198788,Sample_1,A,A,0.0
1:110201112,Sample_1,A,A,0.0042844056
1:110201667,Sample_1,B,B,1.0
1:110202904,Sample_1,A,A,0.0027067661
1:110203240,Sample_1,B,B,1.0
1:110203911,Sample_1,B,B,1.0
1:110206675,Sample_1,A,A,0.0
```

### Procedure

1. Compares samples from `fam` files and intensity files in the directory.
2. Finds markers located in the autosomes.
3. Compute allele frequencies for autosomal markers using _Plink_.
4. Runs [BAFRegress](https://genome.sph.umich.edu/wiki/BAFRegress) using
   autosomal markers.

### Output files

Here is a comprehensive list of all the possible output files for each of the
steps described above.

!!! note
    The output prefix shown is the one by default (_i.e._ `contamination`).

`contamination.to_extract`
: The autosomal markers that will be used by bafRegress.

`contamination.frq`
: The frequency of each of the autosomal markers.

`contamination.bafRegress`
: The bafRegress results for each of the tested sample. See the
[original documentation](https://genome.sph.umich.edu/wiki/BAFRegress#Interpreting_Results)
for results interpretation.

`contamination.contaminated_samples`
: A file containing the list of contaminated samples.
