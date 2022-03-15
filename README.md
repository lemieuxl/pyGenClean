[![PyPI version](https://badge.fury.io/py/pyGenClean.svg)](http://badge.fury.io/py/pyGenClean)


# pyGenClean - Automated Data Clean Up #

`pyGenClean` is an informatics tool to facilitate and standardize the genetic
data clean up pipeline with genotyping array data. In conjunction with a source
batch-queuing system, the tool minimizes data manipulation errors, it
accelerates the completion of the data clean up process and it provides
informative graphics and metrics to guide decision making for statistical
analysis.

If you use `pyGenClean` in you project, please cite the published paper
describing the tool:

> Lemieux Perreault LP, Provost S, Legault MA, Barhdadi A, Dubé MP (2013)
> pyGenClean: efficient tool for genetic data clean up before association testing.
> *Bioinformatics*, **29**(13): 1704-1705
> [DOI:[10.1093/bioinformatics/btt261](http://dx.doi.org/10.1093/bioinformatics/btt261)]



## Documentation

Documentation is available from
[http://lemieuxl.github.io/pyGenClean/](http://lemieuxl.github.io/pyGenClean/).



## Dependencies ##

Here are the dependencies that must be installed before _pyGenClean_:

*   [Python](http://python.org/) (version 2.7)
*   [numpy](http://www.numpy.org/) (version 1.6.2 or latest)
*   [matplotlib](http://matplotlib.org/) (version 1.2.0 or latest)
*   [scipy](http://www.scipy.org/) (version 0.11.0 or latest)
*   [scikit-learn](http://scikit-learn.org/stable/) (version 0.12.1 or latest)
*   [Jinja2](http://jinja.pocoo.org/docs/dev/) (version 2.8 or latest)
*   [PLINK](http://pngu.mgh.harvard.edu/~purcell/plink/) (1.07)



## Installation ##

For Linux users, we recommend installing `pyGenClean` in a Python
[virtualenv](http://pypi.python.org/pypi/virtualenv) (virtual environment).

`pyGenClean` should work on Windows and MacOS, even though it hasn't been fully
tested for full compatibility. It has been tried on Windows XP (32 bits) and
Windows 7 (64 bits, but with a 32 bits Python 2.7 installation) without known
problems.

For a step by step installation on both Linux and Windows operation systems, see
`pyGenClean` documentation, located [here](http://lemieuxl.github.io/pyGenClean/#installation).
