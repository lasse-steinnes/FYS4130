#include "ising_montecarlo.hpp"

// calculating magnetic moments
void IsingMonteCarlo::magnetic_moment2D(){
  /* Code for magnetization for one specific
  state with periodic boundary conditions (2D)
  Calculating total magnetization, by summing over all spins
  for one specific state */
  m_MagneticMoment = 0;
  for (int i = 0; i < m_L*m_L; i++){
    m_MagneticMoment += S2d[i];
  }
}

void IsingMonteCarlo::magnetic_moment1D(){
}
