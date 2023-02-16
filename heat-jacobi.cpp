#include <iostream>
#include <chrono>
#include <cmath>
#include <mpi.h>

#include "a2-helpers.hpp"

using namespace std;

// distrubute to each process
// point to point communication
// allreduce
// gather

int main(int argc, char **argv)
{
    int max_iterations = 15000;
    double epsilon = 1.0e-3;

    // default values for M rows and N columns
    int N = 12;
    int M = 12;

    process_input(argc, argv, N, M, max_iterations, epsilon);

    /* START: the main part of the code that needs to use MPI */

    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    printf("rank %d, size is %d \n", rank, size);

    int i, j;
    double diffnorm, mpidiffnorm;
    int iteration_count = 0;

    // describes the height of the local matrix
    int M_local = (M / size) + 2; //+2 are the buffer zones, one top one bottom.
    int M_local_real = M / size; //height without buffers


    // this doesn't work and I can't figure out why
    // might have been useful if M%size != 0 
    /* 
    if (rank == size - 1){
        M_local = M_local + (M % size);
        M_local_real = M / size + (M % size);
    }  
    */

    // for gatherv
    // number of elements 
    int local_elements = M_local * N;
    int local_elements_real = M_local_real * N;

    // allocate another 2D array
    Mat U_local(M_local, N); // change: use local sizes with MPI, e.g., recalculate M and N
    Mat W(M_local, N);       // change: use local sizes with MPI

    double time_mpi_1 = MPI_Wtime(); // change to MPI_Wait

    // Init + Boundary
    for (i = 0; i < M_local; ++i){
        for (j = 0; j < N; ++j){
            U_local(i, j) = 0.0;
        }
    }

    // adds 100 to the last inner raw
    if (rank == size - 1){
        for (j = 0; j < N; ++j){
            U_local(M_local - 1, j) = 100.0;
        }
    }

    int receive_counts[size];
    int receive_displs[size];

    for (int i = 0; i < size; ++i)
    {
        receive_counts[i] = local_elements_real;
        receive_displs[i] = local_elements_real *i;
    }

    // End init

    iteration_count = 0;
    do
    {

        // communication between overlapping regions 
        // I have used Sendrecv since it was mentioned in the slides that it avoids deadlocks.
        // BUT I have noticed that using sendrecv i get smaller speedups
        // so I have switched back to Send - Recv
        // it seems that never getting a deadlock using sendrecv has an impact on performance

        // raw data measurements in report

        
       /* if (rank < size - 1){
            MPI_Sendrecv(&U_local(M_local_real, 0), N, MPI_DOUBLE, rank + 1, 0, &U_local(M_local_real + 1, 0), N, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        if (rank > 0){
            MPI_Sendrecv(&U_local(1, 0), N, MPI_DOUBLE, rank - 1, 1, &U_local(0, 0), N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        } */
        

        if(rank < size - 1){
             MPI_Send (&U_local(M_local_real,0), N, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD);
        }

        if(rank > 0){
            MPI_Recv (&U_local(0, 0), N, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
        }
        

        if(rank > 0){
             MPI_Send (&U_local(1,0), N, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD);
        }

        if(rank < size - 1){
             MPI_Recv (&U_local(M_local_real+1, 0), N, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD,  MPI_STATUS_IGNORE);
        
        }


        /* Compute new values (but not on boundary) */
        iteration_count++;
        diffnorm = 0.0;

        for (i = 1; i < M_local - 1; ++i){
            for (j = 1; j < N - 1; ++j){
                W(i, j) = (U_local(i, j + 1) + U_local(i, j - 1) + U_local(i + 1, j) + U_local(i - 1, j)) * 0.25;
                diffnorm += (W(i, j) - U_local(i, j)) * (W(i, j) - U_local(i, j));
            }
        }

        // Only transfer the interior points
        for (i = 1; i < M_local - 1; ++i)
            for (j = 1; j < N - 1; ++j)
                U_local(i, j) = W(i, j);


        // computed the diffnorm from all processes
        MPI_Allreduce(&diffnorm, &mpidiffnorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        mpidiffnorm = sqrt(mpidiffnorm);

    } while (epsilon <= mpidiffnorm && iteration_count < max_iterations);


    // calculates elapsed time in the MPI part
    double rank_elapsed_time = MPI_Wtime() - time_mpi_1;
    double rank_time_max = -1;
    MPI_Reduce(&rank_elapsed_time, &rank_time_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    

    // gathers the local matrix from each process into a new "large" matrix
    // ignores the buffers
    Mat U(M_local_real*size, N);
    MPI_Gatherv(&U_local(1, 0), local_elements_real, MPI_DOUBLE, &U(0, 0), receive_counts, receive_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Finalize();

    /* END: the main part of the code that needs to use MPI */

    // verification only if rank == 0
    if (rank == 0){

        // well, the matrix from Gatherv is delayed by one row and I can't figure out
        // how to solve that inside the Gatherv without causing a segmentation fault.
        // this is a workaround but it works.

        Mat U_final(M, N);
        for (int i = 0; i < M_local_real * size; ++i){
            for (int j = 0; j < N; ++j){
                if (i == M_local_real * size - 1){
                    U_final(i, j) = 100;
                }
                else{
                    U_final(i, j) = U(i + 1, j);
                }
            }
        }

        cout << "Computed (MPI) in " << iteration_count << " iterations and " << rank_time_max << " seconds." << endl;

        Mat U_sequential(M, N); // init another matrix for the verification

        // start time measurement for the sequential version
        auto ts_1 = chrono::high_resolution_clock::now();
        heat2d_sequential(max_iterations, epsilon, U_sequential, iteration_count);
        auto ts_2 = chrono::high_resolution_clock::now();

        // print time measurements

        cout << "Computed (sequential) in " << iteration_count << " iterations and " << chrono::duration<double>(ts_2 - ts_1).count() << " seconds." << endl;
        cout << "Verification: " << (verify(U_final, U_sequential) ? "OK" : "NOT OK") << std::endl;

        save_to_disk(U_sequential, "heat2d.txt");
    }

    

    return 0;
}