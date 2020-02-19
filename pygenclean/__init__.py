"""pyGenClean is a tool to facilitate and standardize the genetic QC."""


from .sex_check import sex_check, intensity_plot, baf_lrr_plot


QC_MODULES = [
    ("sex-check", sex_check),
    ("intensity-plot", intensity_plot),
    ("baf-lrr-plot", baf_lrr_plot),
]
