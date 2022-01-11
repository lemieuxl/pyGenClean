"""List of tools available for pyGenClean"""


from .plate_bias import plate_bias
from .related_samples import related_samples
from .sex_check import sex_check, intensity_plot, baf_lrr_plot
from .flag_maf import flag_maf
from .sample_call_rate import sample_call_rate
from .marker_call_rate import marker_call_rate
from .flag_hw import flag_hw
from .hetero_hap import hetero_hap
from .nocall_hetero import nocall_hetero, heterozygosity_plot


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
}


qc_sub_modules = {
    "sex_check": {
        "intensity_plot": intensity_plot,
        "baf_lrr_plot": baf_lrr_plot,
    },
    "nocall_hetero": {
        "heterozygosity_plot": heterozygosity_plot,
    },
}
