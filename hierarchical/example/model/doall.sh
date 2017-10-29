#! /bin/bash -f
ll=loglike.txt
rr=result.txt


python   my_model.py   params.txt data.txt

python   log_like.py    params.txt 4
