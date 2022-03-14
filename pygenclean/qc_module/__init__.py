"""List of tools available for pyGenClean"""


from .plate_bias import plate_bias
from .related_samples import related_samples, merge_related_samples
from .sex_check import (sex_check, intensity_plot, baf_lrr_plot,
                        intensity_viewer)
from .flag_maf import flag_maf
from .sample_call_rate import sample_call_rate
from .marker_call_rate import marker_call_rate
from .flag_hw import flag_hw
from .hetero_hap import hetero_hap
from .nocall_hetero import nocall_hetero, heterozygosity_plot
from .ethnicity import ethnicity, find_outliers, plot_mds, plot_eigenvalues


qc_modules = {
    "plate_bias": plate_bias,
    "related_samples": related_samples,
    "sex_check": sex_check,
    "flag_maf": flag_maf,
    "sample_call_rate": sample_call_rate,
    "marker_call_rate": marker_call_rate,
    "flag_hw": flag_hw,
    "hetero_hap": hetero_hap,
    "nocall_hetero": nocall_hetero,
    "ethnicity": ethnicity,
}


qc_sub_modules = {
    "sex_check": {
        "intensity_plot": intensity_plot,
        "baf_lrr_plot": baf_lrr_plot,
        "intensity_viewer": intensity_viewer,
    },
    "nocall_hetero": {
        "heterozygosity_plot": heterozygosity_plot,
    },
    "related_samples": {
        "merge_related_samples": merge_related_samples,
    },
    "ethnicity": {
        "find_outliers": find_outliers,
        "plot_mds": plot_mds,
        "plot_eigenvalues": plot_eigenvalues,
    }
}
