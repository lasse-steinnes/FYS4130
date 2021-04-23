#include "ising_montecarlo.hpp"
#include <random>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <fstream>

// solver for the isingmodel which gets the expectation values
// time to write the solve algo :)

// solving with monte carlo simulation --> see solve from fys4150
using namespace std;
using namespace chrono;
double * IsingMonteCarlo::solve2D(int calibration, int MC, int rank){
  /*
  Find the expectation values for a 2D steady state Ising Model
  for a given number of cycles
  */

  // calibration: number of calibration cycles
  // returns expectation values
  // calculates all cycles
  // sends in the indices suggested if metropolis gives true
  // update expectation values and flip

  // setting up needed parameters
  m_MC = MC; // number of Monte carlo cycles
  m_rank = rank; // rank of thread in parallelization

  m_calibration = calibration;
  m_L2 = m_L*m_L;
  int n_expv = 3;
  double *exp_values = new double[n_expv]; // The final expectation values
  /*
  Mersenne twister random generator suggest
  flipping of spin with random index. PS: indices are thereafter mapped;
  */

  // seed once in the beginning
  int rd = chrono::high_resolution_clock::now().time_since_epoch().count()+ m_rank; // <--  for parallellization;
  mt19937_64 gen_i(rd);      // seeded with rd
  uniform_int_distribution<> distribution_i(1, (m_L)); // Choose uniform distr. with range 1,(m_L) (unsigned integer)

  int sd = chrono::high_resolution_clock::now().time_since_epoch().count() + m_rank; //  <--  for parallellization;
  mt19937_64 gen_j(sd);     // seeded with sd
  uniform_int_distribution<> distribution_j(1, (m_L)); // Choose uniform distr. with range 1,(m_L)

  //int td = chrono::high_resolution_clock::now().time_since_epoch().count() + m_rank; //<--  for parallellization;
  //m_gen.seed(td);

  double exp_val_M, exp_val_M2;  //expectation values for magnetization and magnetization squared
  double exp_val_Mabs;   //expectation value for mean absolute value of magnetization

  ofstream file_expv; // expectation values file, to close
  open_exp_vals_to_file(file_expv); // opens file to be written

  int rand_i, rand_j;
  for (int temp = 0; temp < m_nT; temp++){
    // Initialize to zero for each temperature
    // state of S at in last MC cycle for previous temp
    exp_val_M = 0.0; exp_val_M2 = 0.0;
    exp_val_Mabs = 0.0;
    m_beta = 1./((double) m_T[temp]);
    m_p = 1.0 - exp(-2*m_beta); // J set to 1, k set to 1

    magnetic_moment2D(); // calculate initial magnetic moment

    // Calibration cycles: run through some percentage of samples before
    // adding energies and magnetic moment to expectation values and variances
    for (int c = 0; c < m_calibration; c++){ // cluster instead of single flip
        rand_i =  distribution_i(gen_i); // Draw index i on physical mesh, suggest flip
        rand_j =  distribution_j(gen_j); // Draw index j on physical mesh, suggest flip
        //cout << "idx :" << rand_i  << " " << rand_j << "\n";
        wolff_sampling2D(rand_i,rand_j);
    }

    // cycles contributing to mean and variance
    for (int c = m_calibration; c < m_MC; c++){ // cluster instead of single flip
        rand_i =  distribution_i(gen_i); // Draw index i on physical mesh, suggest flip
        rand_j =  distribution_j(gen_j); // Draw index j on physical mesh, suggest flip
        wolff_sampling2D(rand_i,rand_j);
      //adding expectation values from each cycle
      exp_val_M += m_MagneticMoment2d;
      exp_val_M2 += m_MagneticMoment2d*m_MagneticMoment2d;
      exp_val_Mabs += fabs(m_MagneticMoment2d);
    }
    //Get final expectation value over all cycles for this temperature: Dividing the sum with number
    //of MC cycles m_MC to get expectation values.

    int c_contrib = m_MC - m_calibration; // cycles contributing to expectation values

    exp_values[0] = exp_val_M/((double) c_contrib);
    exp_values[1] = exp_val_Mabs/((double) c_contrib);

    exp_val_M2 = exp_val_M2/((double) c_contrib);

    //Calculating variance for energy and magnetization
    //Finding specific heat m_Cv and suceptibility m_xi
    double varianceM = exp_val_M2 - exp_values[0]*exp_values[0];
    m_xi = varianceM/((double) m_T[temp]);
    exp_values[2] = m_xi;

    // Scaling by L^2
    for(int i= 0; i < n_expv ;i++){
         exp_values[i] = ((double) 1/m_L2)*exp_values[i];
    }

    // L^4 for variance
    varianceM = ((double) 1/(m_L2*m_L2))*varianceM;

    write_exp_vals_to_file(exp_values,file_expv,temp,varianceM);
  }
  file_expv.close();

  return exp_values; // return something
}
