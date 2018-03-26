This example shows how to couple an external code with TMCMC. In this example the external code is a simple python code but you can substitute it with any code.

1. go to ./source and run make
2. copy 'sample' to the 'example' folder
3. go to the 'example' folder and run './sample'

For each model evaluation needed from TMCMC a folder is created (tmp_dir.\*). The contents of  the folder 'model' are copied into tmp_dir.\* and the doall.sh script is being called. In this  script the model function 'my_model.py' as well as the likelihood function ('log_like.py')  are being called. You can substitute the 'my_model.py' and the 'data.dat' in the 'model' folder  with your favorite model and data file.


At the end of the model evaluation the tmp_dir.\* is not being deleted. If you want to remove
it you have to recompile the code with
```
	#define REMOVEDIRS  1
```
in the source/fitfun.c
