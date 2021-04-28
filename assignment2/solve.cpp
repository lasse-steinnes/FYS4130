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
double * IsingMonteCarlo::solve2D(int calibration, int MC, int N_bins, int rank){
  /*
  Find the expectation values for a 2D steady state Ising Model
  for a given number of cycles
  */

  // calibration: number of calibration cycles
  // returns expectation values
  // calculates all cycles
  // sends in the indices suggested to wolff sampling
  // update expectation values

  // setting up needed parameters
  m_MC = MC; // number of Monte carlo cycles
  m_rank = rank; // rank of thread in parallelization

  m_calibration = calibration;
  m_L2 = m_L*m_L;
  int n_expv = 3;
  double *exp_values = (double *) calloc(n_expv, sizeof(*exp_values));
  /*
  Mersenne twister random generator suggest
  flipping of spin with random index. PS: indices are thereafter mapped;
  */

  // seed once in the beginning
  int rd = chrono::high_resolution_clock::now().time_since_epoch().count()+ m_rank; // <--  for parallellization;
  mt19937_64 gen_i(rd);      // seeded with rd
  uniform_int_distribution<> distribution_i(0, (m_L-1)); // Choose uniform distr. with range 1,(m_L) (unsigned integer)

  int sd = chrono::high_resolution_clock::now().time_since_epoch().count() + m_rank; //  <--  for parallellization;
  mt19937_64 gen_j(sd);     // seeded with sd
  uniform_int_distribution<> distribution_j(0, (m_L-1)); // Choose uniform distr. with range 1,(m_L)

  int td = chrono::high_resolution_clock::now().time_since_epoch().count() + m_rank; //<--  for parallellization;
  m_gen.seed(td);

  double exp_val_M, exp_val_M2;  //expectation values for magnetization and magnetization squared
  double exp_valM4;
  double M2;
  double gamma;

  ofstream file_expv; // expectation values file, to close
  open_exp_vals_to_file(file_expv); // opens file to be written

  int rand_i, rand_j;
  for (int temp = 0; temp < m_nT; temp++){
    // Initialize to zero for each temperature
    // state of S at in last MC cycle for previous temp
    magnetic_moment2D(); // calculate initial magnetic moment
    double T = m_T[temp];
    m_beta = 1./((double) T);
    m_p = 1.0 - exp(-2.0*m_beta); // J set to 1, k set to 1

    // Calibration cycles: run through some percentage of samples before
    // adding magnetic moment to expectation values and variances
    for (int c = 0; c < m_calibration; c++){ // cluster instead of single flip
      // clear the cluster array
      for (int i = 0; i < m_L; i++){
        for (int j = 0; j < m_L; j++){
          cluster[i*m_L + j] = false; // false
          }
        }

        rand_i =  distribution_i(gen_i); // Draw index i on physical mesh, suggest flip
        rand_j =  distribution_j(gen_j); // Draw index j on physical mesh, suggest flip
        growCluster2D(rand_i, rand_j, S2d[rand_i*m_L + rand_j]);
    }

    // cycles contributing to mean and variance
    for (int i_bins = 0; i_bins < N_bins; i_bins++){
      exp_val_M = 0.0; exp_val_M2 = 0.0;
      exp_valM4 = 0.0;

    for (int c = 0; c < m_MC; c++){ // cluster instead of single flip
      // no cluster defined so clear the cluster array
      for (int i = 0; i < m_L; i++){
        for (int j = 0; j < m_L; j++){
          cluster[i*m_L + j] = false; // false
          }
        }

        rand_i =  distribution_i(gen_i); // Draw index i on physical mesh, suggest flip
        rand_j =  distribution_j(gen_j); // Draw index j on physical mesh, suggest flip
        growCluster2D(rand_i, rand_j, S2d[rand_i*m_L + rand_j]);

        //adding expectation values from each cycle
        M2 = m_MagneticMoment*m_MagneticMoment;
        exp_val_M += m_MagneticMoment;
        exp_val_M2 += M2;
        exp_valM4 += M2*M2;

    }
    //Get final expectation value over all cycles for this temperature: Dividing the sum with number
    //of MC cycles m_MC to get expectation values.

    exp_values[0] += ((double) 1/m_MC)*exp_val_M;
    exp_values[1] += ((double) 1/m_MC)*exp_val_M2;
    exp_values[2] +=  ((double) 1/m_MC)*exp_valM4;
    }

    // intrinsic adjust (prop to L) and total mean (divide by bins)
    exp_values[1] = 1./((double) m_L2)*exp_values[1]; // m^2
    exp_values[2] = 1./((double) m_L2*m_L2*m_L2)*exp_values[2]; // m^4

    // Scaling by L^2
    for(int i= 0; i < n_expv ;i++){
         exp_values[i] = 1./((double) m_L2*N_bins)*exp_values[i];
    }

    gamma = exp_values[2]/((double) exp_values[1]*exp_values[1]);
    write_exp_vals_to_file(exp_values,file_expv,T,gamma);
  }
  file_expv.close();
  return exp_values; // return
}

double * IsingMonteCarlo::solve1D(int calibration, int MC, int N_bins, int rank){
  /*
  Find the expectation values for a 1D steady state Ising Model
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
  int n_expv = m_L;
  double *exp_values = new double[n_expv]; // The final expectation values

  /*
  Mersenne twister random generator suggest
  flipping of spin with random index. PS: indices are thereafter mapped;
  */

  // seed once in the beginning
  int rd = chrono::high_resolution_clock::now().time_since_epoch().count()+ m_rank; // <--  for parallellization;
  mt19937_64 gen_i(rd);      // seeded with rd
  uniform_int_distribution<> distribution_i(0, (m_L-1)); // Choose uniform distr. with range 1,(m_L) (unsigned integer)

  int sd = chrono::high_resolution_clock::now().time_since_epoch().count() + m_rank; //  <--  for parallellization;
  mt19937_64 gen_j(sd);     // seeded with sd
  uniform_int_distribution<> distribution_j(0, (m_L-1)); // Choose uniform distr. with range 1,(m_L)

  int td = chrono::high_resolution_clock::now().time_since_epoch().count() + m_rank; //<--  for parallellization;
  m_gen.seed(td);

  int rand_i;
  for (int temp = 0; temp < m_nT; temp++){
    // Initialize to zero for each temperature

    magnetic_moment1D(); // calculate initial magnetic moment
    m_beta = 1./((double) m_T[temp]);
    m_p = 1.0 - exp(-2.0*m_beta); // J set to 1, k set to 1

    cout << "initial magnetic moment: " << m_MagneticMoment<<"\n" ;
    // Calibration cycles: run through some percentage of samples before
    // adding magnetic moment to expectation values and variances
    for (int c = 0; c < m_calibration; c++){ // cluster instead of single flip
      // no cluster defined so clear the cluster array
      for (int i = 0; i < m_L; i++){
          cluster1D[i] = false; // false
        }

        rand_i =  distribution_i(gen_i); // Draw index i on physical mesh, suggest flip
        growCluster1D(rand_i, S1d[rand_i]);
    }

    // cycles contributing to mean and variance
    for (int i_bins = 0; i_bins < N_bins; i_bins++){
      // zero out each time
      for (int i = 0; i < m_L; i++){
          spin_r[i] = 0;
          spin_0r[i] = 0;
        }


    for (int c = 0; c < m_MC; c++){ // cluster instead of single flip
      // no cluster defined so clear the cluster array
      for (int i = 0; i < m_L; i++){
          cluster1D[i] = false; // false
        }

        rand_i =  distribution_i(gen_i); // Draw index i on physical mesh, suggest flip
        growCluster1D(rand_i, S1d[rand_i]);
        get_observables();
    }
    //Get final expectation value over all cycles for this temperature: Dividing the sum with number
    //of MC cycles m_MC to get expectation values.
    get_corr(exp_values);
    }

    // intrinsic adjust (prop to L) and total mean (divide by bins)
    // Scaling by L^2 bc multiplication of terms
    for(int i= 0; i < n_expv ;i++){
         exp_values[i] = 1./((double) m_L2*N_bins)*exp_values[i];
    }
    write_corr(exp_values, temp);
  }
  return exp_values; // return correlation values
}

void IsingMonteCarlo::get_observables(){
  for (int i = 0; i < m_L; i++){
    spin_r[i] += S1d[i];
    spin_0r[i] += S1d[0]*S1d[i];
  }
}

void IsingMonteCarlo::get_corr(double *arr){
  // avg
  for (int i = 0; i < m_L; i++){
    spin_r[i] = ((double) 1/m_MC)*spin_r[i];
    spin_0r[i]= ((double) 1/m_MC)*spin_0r[i];
    arr[i] += spin_0r[i] - spin_r[0]*spin_r[i];
  }
}
