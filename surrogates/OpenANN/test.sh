# 1. HIMMELBLAU
make clean; make dim=2
rm -f myNN.txt
unset THRESHOLD;	./mytest himmelblau_2d_5x1024.txt 1024 0 himmelblau_2d_5x1024.txt 1024 1024
export THRESHOLD=-100;	./mytest himmelblau_2d_5x1024.txt 1024 0 himmelblau_2d_5x1024.txt 1024 1024
