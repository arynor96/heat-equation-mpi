# Solving	Heat	Equation	in	2D

When comparing my MPI implementation to the sequential Heat Equation in 2D, my MPI implementation achieves a 12.39 speedup when using 96 threads.

#### Communication between processes

In P0 the process has to send the last “real” row to the P1 process. At the same time, P0 has to receive from P1, P1’s first “real" row. P0’s first “real” row is not sent above, since there is no process above P0, only below. That’s why all the sends of the first “real” row happen only if the rank is higher than 0, and they have as destination the process below (rank + 1).

In the very last process (Pn-1), the last “real” row can’t be sent anywhere, since there are no processes below. That’s why the last “real” row is sent only if the rank is smaller than size-1.

#### Performance

The combination of MPI_Send and MPI_Recv has a better speedup than MPI_Sendrecv, but I am not sure why. Most likely, it is related to the fact that when using MPI_Sendrecv send and receive happen simultaneously and that means two blocking operations happen simultaneously. That will also explain why after using 4 nodes, the performance is no longer getting better but worse. The amount of simultaneously blockings affects the performance when the number of processes gets even higher.

With the combination of MPI_Send and MPI_Recv, as seen on the graph and raw data, the performance gets better every time the number of processes get higher. With a 2688x4096 matrix and 1000 iterations the speedup is the highest, but after using 4 nodes and 64 the speedup gain is a bit smaller but still far from flattening. When using the same matrix but with 2000 iterations, the speedup gain looks to be more constant between 1 and 6 nodes, around x1.4 for each extra node.

With a 1152x1152 matrix, and 1000 iterations, the speedup is almost the same when using 1,2 or 3 nodes. It gets noticeable better when using 4 nodes. But, out of the all three cases this one has the smallest speedup.

A conclusion might be that the bigger the matrix is, the higher the speedup will be. The number of blocking communications are a bottleneck in my implementation.
