#include "ising_montecarlo.hpp"
#include <random> // To get access to the mersenne twister random generator
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace chrono;

void IsingMonteCarlo::init1D(int L, double T0, double T_end, int n_T){
  // make a string of spin configs
  m_nT = n_T; // number of temperatures to loop over
  m_L = L; // number of spins along a given axis (square 2D system)

  S1d = new int[L];
  cluster1D = new bool[m_L];
  spin_r = new double[m_L];
  spin_0r = new double[m_L];

  if(T0 >= 1.1) {        //Temperature check
    for(int i = 0; i < L; i++) {    //If the temperature is greater than 1.1,
      if(m_check < 0.5) {             //the lattice is filled with random spins.
        S1d[i] = -1;
      }else {
        S1d[i] = 1;
      }
      draw_acceptance();
    }

  }else {
    for(int i = 0; i < L; i++) {    //If the temperature is smaller than 1.1,
      if(m_check < 0.5) {             //the lattice is filled with either only
        S1d[i] = -1;                    //positive spins, or only negative.
      }else {
        S1d[i] = 1;
      }
    }
  }


  // filling in the T vector to be used in MC
  m_T = new double[m_nT];
  double dT;
  if (n_T > 1){
    dT = (T_end-T0)/((double) (n_T-1));
  }else{
    dT = 0;
  }

  for (int i = 0; i < m_nT; i++){
    m_T[i] = T0 + i*dT;
  }

};


void IsingMonteCarlo::init2D(int L, double T0, double T_end, int n_T){
  // make a square spin lattice with
  // ghost cells to enforce periodic boundary conditions
  // map is index mapping vector of length L + 2

  // allocate 2-D array for spin cluster labels

  m_nT = n_T; // number of temperatures to loop over
  m_L = L; // number of spins along a given axis (square 2D system)
  cluster = new bool[m_L*m_L];

  // Random spin configuration
  //Setting up lattice of L*L elements
  S2d = new int[L*L];
  draw_acceptance();    //Getting random number
  if(T0 >= 1.1) {        //Temperature check
    for(int i = 0; i < L*L; i++) {    //If the temperature is greater than 1.1,
      if(m_check < 0.5) {             //the lattice is filled with random spins.
        S2d[i] = -1;
      }else {
        S2d[i] = 1;
      }
      draw_acceptance();
    }

  }else {
    for(int i = 0; i < L*L; i++) {    //If the temperature is smaller than 1.1,
      if(m_check < 0.5) {             //the lattice is filled with either only
        S2d[i] = -1;                    //positive spins, or only negative.
      }else {
        S2d[i] = 1;
      }
    }
  }

  m_T = new double[m_nT];
  double dT;
  if (n_T > 1){
    dT = (T_end-T0)/((double) (n_T-1));
  }else{
    dT = 0;
  }

  for (int i = 0; i < m_nT; i++){
    m_T[i] = T0 + i*dT;
  }
}
