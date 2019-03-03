
This is a package of procedures for fitting models to infrared spectral energy distributions (SEDs) of galaxies. 

The files makeplots2_rieke.py and makeplots2_chary.py are scripts that can be run to read galaxy templates, read theoretical models from the Draine & Li (2007) library, construct Monte Carlo realizations of the templates with errors, fit models to the data, and quantify the biases associated with fitting one model at a time vs fitting linear combinations of models.

To get these scripts to run, you will need to:

- download the Draine & Li 2007 model tarfile (not included here since it is large) and untar it in the draine_li_2007 subdirectory. See the Readme.models file. See composite_draineli_sed.py, which constructs composite models.
- fix some of the filenames in the read*.py files, which are currently pointing at different directories on my computer. They should be reading from the data, filters, and draine_li_2007 subdirectories (to be fixed)

- make subdirectories for the scripts to save output plots it (see "plotdir" in the scripts)

This work is currently unpublished but will appear in a conference proceeding; please contact me for how to reference it. 

Benjamin Weiner, bjw@as.arizona.edu
