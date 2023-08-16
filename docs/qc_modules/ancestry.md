# Ancestry

The _ancestry_ QC module consists of a main script and three companion scripts.
It performs a PCA analysis (with or without reference populations) and generates
multiple ancestry plots.

Use the following command to access the multiple scripts of the _ancestry_ QC
module.

```shell-session
$ pyGenClean ancestry --help
usage: pyGenClean ancestry [-h] [-v]
                           {run,find-outliers,plot-mds,plot-eigenvalues} ...

Checks sample's ancestry using reference populations.

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

subcommands:
  Below is a list of tools in the ancestry module. Note that 'run' executes
  the main ancestry pipeline.

  {run,find-outliers,plot-mds,plot-eigenvalues}
    run                 Checks sample's ancestry using reference populations.
    find-outliers       Finds outliers in SOURCE from CEU samples.
    plot-mds            Create an MDS plot.
    plot-eigenvalues    Plot eigenvalues (scree plot).
```

## Main script

The main script (access by using the `run` subcommand) performs a PCA analysis
of the source dataset (with or without reference populations).

```shell-session
$ pyGenClean ancestry run --help
usage: pyGenClean ancestry run [-h] --bfile FILE [--skip-ref-pops]
                               [--ceu-bfile FILE] [--yri-bfile FILE]
                               [--jpt-chb-bfile FILE] [--plink-1.07]
                               [--nb-threads N] [--min-nb-snp INT]
                               [--indep-pairwise STR STR STR] [--maf FLOAT]
                               [--nb-components INT]
                               [--find-outliers-format FORMAT]
                               [--outliers-of POP] [--multiplier FLOAT]
                               [--find-outliers-xaxis COMPONENT]
                               [--find-outliers-yaxis COMPONENT]
                               [--plot-mds-format FORMAT]
                               [--plot-mds-title STRING]
                               [--plot-mds-xaxis STRING]
                               [--plot-mds-yaxis STRING] [--create-scree-plot]
                               [--plot-eigen-title TITLE] [--out FILE]

Checks sample's ancestry using reference populations.

options:
  -h, --help            show this help message and exit

Input File:
  --bfile FILE          The input file prefix (will find the plink binary
                        files by appending the prefix to the .bim, .bed and
                        .fam files, respectively.
  --skip-ref-pops       Perform the MDS computation, but skip the three
                        reference panels.
  --ceu-bfile FILE      The input file prefix (will find the plink binary
                        files by appending the prefix to the .bim, .bed and
                        .fam files, respectively) for the CEU population
  --yri-bfile FILE      The input file prefix (will find the plink binary
                        files by appending the prefix to the .bim, .bed and
                        .fam files, respectively) for the YRI population
  --jpt-chb-bfile FILE  The input file prefix (will find the plink binary
                        files by appending the prefix to the .bim, .bed and
                        .fam files, respectively) for the JPT-CHB population

Options:
  --plink-1.07          Use original Plink (version 1.07)
  --nb-threads N        The number of threads for this analysis (no effect
                        when using plink 1.07). [1]
  --min-nb-snp INT      The minimum number of markers needed to compute IBS
                        values. [Default: 8000]
  --indep-pairwise STR STR STR
                        Three numbers: window size, window shift and the r2
                        threshold. [default: ['50', '5', '0.1']]
  --maf FLOAT           Restrict to SNPs with MAF >= threshold. [default:
                        0.05]
  --nb-components INT   The number of component to compute. [default: 10]

Outlier Options:
  --find-outliers-format FORMAT
                        The output file format (png, ps, or pdf formats are
                        available). [default: png]
  --outliers-of POP     Finds the outliers of this population. [default: CEU]
  --multiplier FLOAT    To find the outliers, we look for more than x times
                        the cluster standard deviation. [default: 1.9]
  --find-outliers-xaxis COMPONENT
                        The component to use for the X axis. [default: C1]
  --find-outliers-yaxis COMPONENT
                        The component to use for the Y axis. [default: C2]

MDS Plot Options:
  --plot-mds-format FORMAT
                        The output file format (png, ps, pdf, or X11 formats
                        are available). [default: png]
  --plot-mds-title STRING
                        The title of the MDS plot. [default: C2 in function of
                        C1 - MDS]
  --plot-mds-xaxis STRING
                        The component to use for the X axis. [default: C1]
  --plot-mds-yaxis STRING
                        The component to use for the Y axis. [default: C2]

Scree Plot Options:
  --create-scree-plot   Computes Eigenvalues and creates a scree plot.
  --plot-eigen-title TITLE
                        The main title of the scree plot [EIGENSOFT results]

Output File:
  --out FILE            The prefix of the output files. [default: ancestry]
```

### Input files

This module uses PLINK’s binary file format (`bed`, `bim` and `fam` files) for
the source data set (the data of interest) and, optionally, three sets of binary
files for the reference panels (_CEU_, _YRI_ and _JPG-CHB_).

### Procedure

1. Finds the overlapping markers between the three reference panels and
   the source panel.
2. Extract the required markers from all the data sets.
3. Renames the reference panel's marker names to that they are the same as
   the source panel (for all populations).
4. Combines the three reference panels together.
5. Compute the frequency of all the markers from the reference and the
   source panels.
6. Finds the markers to flip in the reference panel (when compared to the
   source panel).
7. Excludes the unflippable markers from the reference and the source
   panels.
8. Flips the markers that need flipping in their reference panel.
9. Combines the reference and the source panels.
10. Runs part of the related samples QC module on the combined data set to
    generate the `genome` file.
11. Creates the `mds` file from the combined data set and the result of
    previous step.
12. Creates the population file.
13. Plots the `mds` values.
14. Finds the outliers of a given reference population's cluster.
15. If required, computes the Eigenvalues using `smartpca`.
16. If required, creates a scree plot from `smartpca` resutls.

### Output files

The output files of each of the steps described above are as follow (note that
the output prefix shown is the one by default, _i.e._ ancestry):

`ancestry.ref_snp_to_extract`

: The list of markers to extract from the reference panels.

`ancestry.source_snp_to_extract`

: The list of markers to extract from the source panel.

`ancestry.update_names`

: The updated names of the marker in the reference panels, so that they match
  with the names in the source panel.

`ancestry.reference_panel.CEU`

: The data set containing the extracted markers from the CEU reference
  population.

`ancestry.reference_panel.YRI`

: The data set containing the extracted markers from the YRI reference
  population.

`ancestry.reference_panel.JPT-CHB`

: The data set containing the extracted markers from the JPG-CHB reference
  population.

`ancestry.source_panel.ALL`

: The data set containing the extracted markers from the source population.

`ancestry.reference_panel.ALL.files_to_merge`

: The file required by Plink to merge more than two data sets together.

`ancestry.reference_panel.ALL`

: The data set containing the merged data sets of the three reference
  population.

`ancestry.reference_panel.ALL.rename`

: The data set after markers have been renamed in the reference panels.

`ancestry.reference_panel.ALL.rename.frequency`

: The frequencies of the markers in the reference panels.

`ancestry.source_panel.ALL.frequency`

: The frequencies of the markers in the source panels.

`ancestry.snp_to_flip_in_reference`

: The list of markers to flip in the reference panels.

`ancestry.snp_to_remove`

: The list of markers to remove because they are not comparable to the markers
  in the source panel, even after trying to flip them.

`ancestry.reference_panel.ALL.rename.cleaned`

: The data set after the markers found in the previous step are excluded from
  the reference panels.

`ancestry.source_panel.ALL.cleaned`

: The data set after the markers found in the previous step are excluded from
  the source panel.

`ancestry.reference_panel.ALL.rename.cleaned.flipped`

: The data set after markers from the reference panels were flipped so that they
  become comparable with the source panel.

`ancestry.final_dataset_for_genome.files_to_merge`

: The file required by Plink to merge more than two data sets together.

`ancestry.final_dataset_for_genome`

: The data set containing the merged reference and source panels.

`ancestry.ibs`

: For more information about those files (see the
  [related samples](related_samples.md) QC module for more information).

`ancestry.mds`

: Files containing the MDS values.

`ancestry.population_file`

: The population file required for MDS value plotting.

`ancestry.mds.png`

: The plot of the MDS values (see Figure Initial MDS plot).

`ancestry.before.png`

: The MDS values before outliers detection (see Figure MDS plot before outlier
  detection).

`ancestry.after.png`

: The MDS values after outliers detection for each of the three reference
  populations. The shaded points are the outliers (see Figure MDS plot after
  outlier detection).

`ancestry.outliers.png`

: The MDS values after outliers detection for the selected reference population
  (default is CEU) (see Figure Ethnic outliers).

`ancestry.outliers`

: The list of outliers (excluding the reference populations).

`ancestry.population_file_outliers`

: A population file containing the outliers (to help creating a new MDS plot
  using pyGenClean.PlinkUtils.plot_MDS_standalone).

### Figures

Multiple figures are generated.
