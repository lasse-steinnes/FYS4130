#include "ising_montecarlo.hpp"
#include <chrono>
#include <random>

using namespace std;
using namespace chrono;

void IsingMonteCarlo::draw_acceptance(){ // used to initialize
  /* Mersenne twister random generator suggest
  random number between [0,1) */
  int td = chrono::high_resolution_clock::now().time_since_epoch().count(); // Used to obtain seed
  mt19937_64 gen(td);                                                       // seeded with sd
  uniform_real_distribution<double> distribution(0.0,1.0);                  // creates [0,1)
  m_check =  distribution(gen);                                           // draw acceptance criteria
}; // gen is the argument to distribution
