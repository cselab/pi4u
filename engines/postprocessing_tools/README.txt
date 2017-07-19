#-------------------------- GENERAL SETUP ------------------------------
An R script for plotting 2D projections of samples and estimate the 2D and marginal density is given in the file kdeplotter.R. To make it work:
	o	provide full paths to the curgen_db file and to the helper files of kdeplotter (see comments in the script)
	o	set other user-defined variables (like labels for the plot)
	o	run (with `Rscript kdeplotter.R`)
	o	the picture will appear in the same folder as the curgen_db file and will have the same name with the corresponding extension

#-------------------------- LIBRARIES ----------------------------------
If the "sp" and "MASS" libraries are not installed, run
install.packages("sp")
install.packages("MASS")

#---------------------------- EXAMPLE -----------------------------------
o	run `Rscript kdeplotter.R` from the current directory
o	the output file curgen_db_007.png should be produced
