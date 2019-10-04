
#ifndef MODULES_TEST_TASKS_TEST_MPI_OPS_MPI_H_
#define MODULES_TEST_TASKS_TEST_MPI_OPS_MPI_H_

#include <vector>

double getParallelIntegration(std::vector<double> vec, double* F, double a, double b, int N);
double getSequentialIntegration(std::vector<double> vec, double* F, double l_h);

#endif  // MODULES_TEST_TASKS_TEST_MPI_OPS_MPI_H_