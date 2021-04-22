#ifndef ISINGMONTECARLO_HPP
#define ISINGMONTECARLO_HPP

#include <iostream>
#include <chrono>
#include <random> // To get access to the mersenne twister random generator

using namespace std;

class IsingMonteCarlo{
protected:  // allocate pointers inside main which is accessable to all functions
int m_rank; // rank
int m_nT;   // number of temps to loop over
double m_beta; // a given temperature
int m_L;  // number of spins in 1 dimension
int m_MC; // monte carlo cycles
int *m_map; // mapping
int *S2d; // spin matrix
int *S1d; // spin vector
double m_check; // random number to initialize
int m_MagneticMoment; // total magnetization  (avg  = M/L^2)
double m_p; // probability of drawing one spin (1-exp(-2*m_beta))
// Use this to gen_accept
 /*random number generator, seeded for each rank once */

public:
  void init1D(int L, double T0, int n_T, int MC, int rank); // initiates the 1D spin system
  void init2D(int L, double T0, int n_T, int MC, int rank, int *map); // initiates the 2D spin system
  void magnetic_moment1D();
  void magnetic_moment2D();
  void wolff_sampling2D(int node_i, int node_j); // sampling algo
  void wolff_sampling1D(int node_id);
  void write_spin_to_file(); //writes spin config to file;
  void write_exp_vals_to_file(); //write expectation values to file
  void draw_acceptance(); // Mersenne twister
  void solve(); //uses other functions to solve and get observables
};

#endif
