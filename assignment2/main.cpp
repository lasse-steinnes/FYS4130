// main for assignment 2 in statistical mechanics about the Ising Model
//#define CATCH_CONFIG_RUNNER // This tells Catch to not provide a main()
//#include "catch.hpp"
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
  int n_T, numthreads, dims;
  bool save_over_cycles = false;
  bool save_spin = false;
  int calib;
  int rank = 0 ;

  cout << "Enter integer number of spin particles for each axis:" << " ";
  cin >> L;
  cout << "Enter start point temperature:"  << " ";
  cin >> T_start;
  cout << "Enter an endpoint temperature:"  << " ";
  cin >> T_end;
  cout << "Enter integer number of temperature points to be evaluated:"  << " ";
  cin >> n_T;
  cout << "Enter integer number of MC cycles:"  << " ";
  cin >> MC;
  cout << "Enter integer number of calibration cycles:"  << " ";
  cin >> calib; // eg. 20 000
  cout << "Enter integer number of threads:"  << " ";
  cin >> numthreads;
  cout << "Enter dimension (1/2) ";
  cin >> dims;

  // make a temperature pointer to use
  //double *T = new double[n_T];
  //double dT = (T_end-T_start)/(n_T-1);
  //for (int i = 0; i < n_T; i++){
  //  T[i] = T_start + i*dT;
  //}

  // setup mapping (1D and 2D) -- not sure if this is needed
  int * map = new int[L + 2]; // allocate
  map[0] = L - 1; // first index in map --> last index in physical mesh
  map[L+1] = 0; // last index in map --> first index in physical mesh
  int N_bins = 10;

  if (dims == 1){
    // intialize 1D case
    //Solver.init1D();
  }

  if (dims == 2){
    IsingMonteCarlo Solver;
    // initialize 2D case

  for (int i = 0; i < L; i++){
    map[i+1] = i;
  }

    double *ana = new double[3];

    Solver.init2D(L, T_start,T_end, n_T, map);
    ana = Solver.solve2D(calib, MC, N_bins, rank);

    /*
  int calibration = 20000;
  int L = 2;
  int MC = 1e7;
  double T = 1.0;
  double T2 = 1.0;
  int n = 1;
  double B = 1./T;
  int rank = 0;
    */

    double T = 1.0;
    double B = 1./T;
    double z = 12. + 4*cosh(8*B);
    double mean_abs_M = (8./z)*(exp(8*B) + 2);
    double xi = ((double) 32/(T*z))*(exp(8*B) + 1);

    for (int i = 0; i < 3; i++){
      cout << ana[i] << "\n";
    }

    cout << xi/((double) (2*2*2*2)) << "\n"; // prop var dep L^4
    cout << mean_abs_M/((double) (2*2)) << "\n"; // prop L^2

  }

  if (dims == 0||dims > 2){ // || means or
    cout << "\n\n";
    cout << "Wrong usage, choose correct task options!\n";
    return 1;
    }

  return 0;

}
