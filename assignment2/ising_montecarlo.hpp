#ifndef ISINGMONTECARLO_HPP
#define ISINGMONTECARLO_HPP

#include <iostream>
#include <chrono>
#include <random> // To get access to the mersenne twister random generator
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

class IsingMonteCarlo{
protected:  // allocate pointers inside main which is accessable to all functions
int m_rank; // rank
int m_nT;   // number of temps to loop over
double m_beta; // a given temperature
int m_L, m_L2;  // number of spins in 1 dimension
int m_MC; // monte carlo cycles
double m_xi; //magnetic susceptibility
int *m_map; // mapping
int *S2d; // spin matrix
int *S1d; // spin vector
double *m_T; // temperature vector
double m_check; // random number to initialize
int m_MagneticMoment2d; // total magnetization  (avg  = M/L^2)
int m_MagneticMoment1d; // total magnetization  (avg  = M/L^2)
double m_p; // probability of drawing one spin (1-exp(-2*m_beta))
int m_calibration; // number of calibration cycles

public:
  void init1D(int L, double T0, double T_end, int n_T, int *map); // initiates the 1D spin system
  void init2D(int L, double T0, double T_end, int n_T, int *map); // initiates the 2D spin system
  void magnetic_moment1D();
  void magnetic_moment2D();
  void wolff_sampling2D(int node_i, int node_j); // sampling algo
  void wolff_sampling1D(int node_id);
  void open_exp_vals_to_file(ofstream&file); //writes spin config to file;
  void write_exp_vals_to_file(double *expval,ofstream&file, int temp, double varM); //write expectation values to file
  void draw_acceptance(); // Mersenne twister
  double * solve2D(int calibration, int MC, int rank); //uses other functions to solve and get observables
};

#endif
