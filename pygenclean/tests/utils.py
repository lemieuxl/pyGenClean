"""Utility functions for testing."""


import random
from pathlib import Path
from typing import Optional

import numpy as np
from pyplink import PyPlink


def generate_plink_files(bfile: str, genotypes: np.array,
                         chromosomes: Optional[np.array] = None) -> Path:
    """Generate Plink files randomly."""
    nb_samples = genotypes.shape[1]

    # Creating the FAM file
    with open(bfile + ".fam", "w") as fam:
        for i in range(nb_samples):
            sample_id = f"s{i}"
            sex = random.randint(1, 2)
            print(sample_id, sample_id, 0, 0, sex, -9, file=fam)

    # Generating chromosomes, of none provided (autosome only)
    if chromosomes is None:
        chromosomes = np.random.randint(1, 23, size=genotypes.shape[0])
    chromosomes = np.sort(chromosomes)

    # Creating the BED and BIM files
    with PyPlink(bfile, "w") as bed, open(bfile + ".bim", "w") as bim:
        for i, (chrom, geno) in enumerate(zip(chromosomes, genotypes)):
            bed.write_genotypes(geno)
            print(chrom, f"var{i}", 0, i + 1, "A", "B", sep="\t", file=bim)


def generate_genotypes(nb_samples: int, nb_variants: int) -> np.array:
    """Generate random genotypes."""
    # Generating genotypes
    genotypes = np.empty((nb_variants, nb_samples), dtype=np.int8)
    for i in range(nb_variants):
        genotypes[i] = np.random.binomial(
            n=2, p=random.random() / 2, size=nb_samples,
        )

    # According to the number of samples, MAF could be above 0.5
    to_fix = (np.mean(genotypes, axis=1) / 2) > 0.5
    genotypes[to_fix] = 2 - genotypes[to_fix]

    return genotypes
