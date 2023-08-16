# Welcome to pyGenClean’s documentation

## Introduction

Genetic association studies making use of high throughput genotyping arrays need
to process large amounts of data in the order of millions of markers per
experiment. The first step of any analysis with genotyping arrays is typically
the conduct of a thorough data clean up and quality control to remove poor
quality genotypes and generate metrics to inform and select individuals for
downstream statistical analysis.

`pyGenClean` is a bioinformatics tool to facilitate and standardize the genetic
data clean up pipeline with genotyping array data. In conjunction with a source
batch-queuing system, the tool minimizes data manipulation errors, it
accelerates the completion of the data clean up process and it provides
informative graphics and metrics to guide decision making for statistical
analysis.

`pyGenClean` is a command tool working on both Linux and Windows operating
systems. Its usage is shown below:

```shell-session
$ pyGenClean --help
usage: pyGenClean [-h] [-v] [--debug]
                  {pipeline,plate-bias,related-samples,sex-check,flag-maf,sample-call-rate,marker-call-rate,flag-hw,hetero-hap,nocall-hetero,ancestry,subset,contamination}
                  ...

Runs pyGenClean.

positional arguments:
  {pipeline,plate-bias,related-samples,sex-check,flag-maf,sample-call-rate,marker-call-rate,flag-hw,hetero-hap,nocall-hetero,ancestry,subset,contamination}
    pipeline            The main pyGenClean pipeline.
    plate-bias          Check for plate bias.
    related-samples     Finds related samples according to IBS values.
    sex-check           Check sample's sex using Plink.
    flag-maf            Flag markers with MAF of 0 (monomorphic).
    sample-call-rate    Remove samples with poor call rate.
    marker-call-rate    Remove markers with poor call rate.
    flag-hw             Flags markers failing Hardy Weinberg equilibrium.
    hetero-hap          Remove heterozygous haploid (males on chromosome X).
    nocall-hetero       Clean markers with no call or heterozygous only.
    ancestry            Checks sample's ancestry using reference populations.
    subset              Subset Plink datasets (samples or markers).
    contamination       Check for sample contamination using BAFRegress.

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  --debug               Enable debug logging.
```

The tool consists of multiple standalone scripts that are linked together via a
main script (`pyGenClean`) and a configuration file (the `--conf` option), the
latter facilitating user customization.

The Data clean up protocol schema shows the proposed data cleanup pipeline. Each
box represents a customizable standalone script with a quick description of its
function. Optional manual checks for go-no-go decisions are indicated.

## Installation

`pyGenClean` is a Python package that works on both Linux and Windows operating
systems. It requires a set of Python dependencies and PLINK. Complete
installation procedures are available for both Linux and Windows in the
following sections.

- Linux Installation
- Windows Installation

## Input files

To use `pyGenClean`, two sets of files are required: a set of genotype files and
a configuration file (using the `TOML` format).

### Genotype files

The input files (`--bfile`) of the main program (`pyGenClean`) are _PLINK_'s
binary pedfile format consists of three files with the following extensions:
_BED_, _BIM_ and _FAM_.

For more information about these file formats, have a look at PLINK’s website,
in the Basic usage/data formats section
([https://zzz.bwh.harvard.edu/plink/data.shtml](https://zzz.bwh.harvard.edu/plink/data.shtml)).

### Configuration file

To customized `pyGenClean`, a basic configuration file is required. It tells
which script to use in a specific order. It also sets the different options and
input files, so that the analysis is easy to replicate or modify.

The configuration file uses the [`TOML`](https://toml.io/) format. It consists
of a "_table_" of steps followed by its customization. Lines preceded by a `#`
are comments and are not read by `pyGenClean`.

The following example first removes samples with a missing rate of 10% and more
(section starting at line 1), then removes markers with a missing rate of 2% and
more (section starting at line 6). Finally, it removes the samples with a
missing rate of 2% and more (section starting at line 3).

```toml linenums="1"
[steps.1]
# Removes sample with a missing rate higher than 10%.
module = sample-call-rate
mind   = 0.1

[steps.2]
# Removes markers with a missing rate higher than 2%.
script = marker-call-rate
geno   = 0.02

[steps.3]
# Removes sample with a missing rate higher than 2%.
script = sample-call-rate
mind   = 0.02
```

For a more thorough example, complete configuration files are available for
download at [https://statgen.org](https://statgen.org) and are explained in the
Configuration Files section. For a list of available QC modules, refer to the
List of Modules and their Options.
