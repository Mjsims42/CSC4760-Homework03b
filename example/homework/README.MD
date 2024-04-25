00 command to run:
mpic++ init.cpp -o host
mpirun -n (NodeSize) --oversubscribe host (P) (Q)

01 command to run:
mpic++ init.cpp -o host
mpirun -n (NodeSize) --oversubscribe host (P) (Q) (M)

02 command to run:
make
mpirun -n (NodeSize) --oversubscribe prog.host (M) (N) (P) (Q) (Iterations)

03 command to run:
mpic++ init.cpp -o host
mpirun -n (NodeSize) --oversubscribe host (P) (Q) (M)
