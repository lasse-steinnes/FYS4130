// main for assignment 2 in statistical mechanics about the Ising Model
#include "ising_montecarlo.hpp"
#include <iostream>
#include <omp.h>
#include <stdio.h>
#include <chrono>

using namespace std;
using namespace chrono;

int main(int argc, char const *argv[]){
  int L; int MC;
  double T_start, T_end;
  int num_threads;
  int n_T, dims;
  int calib;
  int rank = 0 ;

  cout << "Enter integer number of spin particles for each axis: ";
  cin >> L;
  cout << "Enter start point temperature: ";
  cin >> T_start;
  cout << "Enter an endpoint temperature: ";
  cin >> T_end;
  cout << "Enter integer number of temperature points to be evaluated: ";
  cin >> n_T;
  cout << "Enter integer number of MC cycles: ";
  cin >> MC;
  cout << "Enter integer number of calibration cycles: ";
  cin >> calib; // eg. 20 000
  cout << "Enter integer number of threads: ";
  cin >> num_threads;
  cout << "Enter dimension (1/2) ";
  cin >> dims;

  // setup mapping (1D and 2D) -- not sure if this is needed
  int N_bins = 10;

  IsingMonteCarlo Solver;
  double *ana = new double[4];

  if (dims == 1){
    // intialize and solve 1D case

    Solver.init1D(L, T_start,T_end, n_T);
    ana = Solver.solve1D(calib, MC, N_bins, rank);

    cout << "numerical: \n";
    for (int i = 0; i < L; i++){
      cout << ana[i] << "\n";
    }
    cout << "Above: Correlation values\n";

  }

  if (dims == 2){
    // initialize 2D case
    // parallellize this

    // define object types
   int temps_i;
   double dT;
   double start, end;
   double *T_vec;
   T_vec = new double[num_threads + 1];
   cout << "num threads: " << num_threads << "\n";

   if (num_threads > 1){
     dT = (T_end-T_start)/((double) (num_threads)); // if num threads 4 --> 5 points
   }else{
     dT = 0;
   }

   for (int i = 0; i < num_threads+1; i++){
     T_vec[i] = T_start + i*dT;
   }


  start = omp_get_wtime();
  omp_set_num_threads(num_threads);
  #pragma omp parallel for private(temps_i)
  for (temps_i = 0; temps_i < num_threads; temps_i++){
    IsingMonteCarlo Solver; // initate class object;
    Solver.init2D(L, T_vec[temps_i],T_vec[temps_i+1], n_T);
    ana = Solver.solve2D(calib, MC, N_bins,  omp_get_thread_num());
    printf("Thread rank: %d\n", omp_get_thread_num());
  }
  end = omp_get_wtime();
  printf("Work took %f seconds\n", end - start);
  }

  if (dims == 0||dims > 2){ // || means or
    cout << "\n\n";
    cout << "Wrong usage, choose correct task options!\n";
    return 1;
    }

  return 0;

}
