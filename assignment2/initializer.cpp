#include "ising_montecarlo.hpp"
#include <random> // To get access to the mersenne twister random generator
#include <iostream>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace chrono;

void IsingMonteCarlo::init1D(int L, double T0, int n_T, int MC, int rank){
  // make a string of spin configs
  m_rank = rank; // rank of thread in parallelization
  m_nT = n_T; // number of temperatures to loop over
  m_L = L; // number of spins along a given axis (square 2D system)
  m_MC = MC; // number of Monte carlo cycles

  S1d = new int[L];

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


};


void IsingMonteCarlo::init2D(int L, double T0, int n_T, int MC, int rank, int *map){
  // make a square spin lattice with
  // ghost cells to enforce periodic boundary conditions
  // map is index mapping vector of length L + 2

  m_map = map; // index mapping
  m_rank = rank; // rank of thread in parallelization
  m_nT = n_T; // number of temperatures to loop over
  m_L = L; // number of spins along a given axis (square 2D system)
  m_MC = MC; // number of Monte carlo cycles

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

}
