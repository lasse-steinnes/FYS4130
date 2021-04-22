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
}; // gen is the argument to dist_accept

/*
// seed once in the beginning
 int rd = chrono::high_resolution_clock::now().time_since_epoch().count()+ m_rank; // <--  for parallellization;
 mt19937_64 gen_i(rd);      // seeded with rd
 uniform_int_distribution<> distribution_i(1, (m_L)); // Choose uniform distr. with range 1,(m_L) (unsigned integer)

 int sd = chrono::high_resolution_clock::now().time_since_epoch().count() + m_rank; //  <--  for parallellization;
 mt19937_64 gen_j(sd);     // seeded with sd
 uniform_int_distribution<> distribution_j(1, (m_L)); // Choose uniform distr. with range 1,(m_L)
 */ // generate the random key node
