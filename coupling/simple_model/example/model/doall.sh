#! /bin/bash -f
ll=loglike.txt
rr=result.txt


module load python

python   my_model.py   params.txt  data.dat

python   log_like.py    params.txt 4
