rm -f curgen_db*.txt;
#mpirun -n 2 ./abc_subsim_torc.exe 2000 curgen_db 0 10 0.05 0 0 2
 mpirun -n 1 ./abc_subsim_torc.exe 2000 curgen_db 0 30 0.05 0 0 2
./display_gen.exe curgen_db 0 2 
./display_gen.exe curgen_db 5 2 
