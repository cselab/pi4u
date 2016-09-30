#-------------------------- GENERAL SETUP ------------------------------
An R script for plotting 2D projections of samples and estimate the 2D and marginal density is given in the file kdeplotter.R. To make it work:
	- provide full paths to the curgen_db file and to the helper files of kdeplotter (see comments in the script)
	- set other user-defined variables (like labels for the plot)
	- run (with source("kdeplotter.R"))
	- the picture will appear in the same folder as the curgen_db file and will have the same name with the corresponding extension

#-------------------------- LIBRARIES ----------------------------------
If the "sp" and "MASS" libraries are not installed, run
install.packages("sp")
install.packages("MASS")

#---------------------------- EXAMPLE -----------------------------------
run source("kdeplotter.R") in the R command line from the current directory
the output file curgen_db_007.png should be produced
