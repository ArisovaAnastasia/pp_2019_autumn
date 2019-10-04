
#include <mpi.h>
#include <vector>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_1/arisova_a_rectangle_integration/rectangle_integration.h"


double getSequentialIntegration(std::vector<double> vec, double* F, double l_h) {
    double t_sum=0;
    for(int i = 0; i < N; i++){
       t_sum += l_h/2 * (F(vec[i]) + F(vec[i+1]));
    }
    return t_sum;
}

double getParallelIntegration(std::vector<double> vec, double* F, double a, double b, int N) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    std::vector<double>global_vector(N + 1);
    double h= (b-a) / (double) N;
    for(int i = 0; i < N; i++){
        global_vector[i] = a + i * h;
    }
    global_vector[4] = b;

     int l_N = (N + 1) /size;

    if (rank == 0) {
        for (int proc = 1; proc < size; proc++) {
            MPI_Send(&global_vector[0] + proc * l_N, l_N,
                        MPI_INT, proc, 0, MPI_COMM_WORLD);
        }
    }

    std::vector<double> local_vector(l_N);
    if (rank == 0) {
        local_vec = std::vector<int>(global_vector.begin(),
                                     global_vector.begin() + l_N);
    } else {
        MPI_Status status;
        MPI_Recv(&local_vector[0], l_N, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }

    int g_sum = 0;
    int l_sum = getSequentialIntegration(local_vector, F, h);
    MPI_Op op_code;
    op_code = MPI_SUM;     
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, op_code, 0, MPI_COMM_WORLD);
    return global_sum;
}
