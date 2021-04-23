// Programmes to write expectation values to file
#include "ising_montecarlo.hpp"


void IsingMonteCarlo::open_exp_vals_to_file(ofstream&file){ // write expectation values after all cycles
  string filename("./Results/Data/expvalues" + to_string(m_MC) + \
                  "-" + to_string(m_L) + "by" + to_string(m_L) + "rank" + to_string(m_rank) + ".txt");
  file.open(filename);
  file    << "T" << setw(25) << "MC_cycles-ac" << setw(25) << "N_spins" << setw(25)\
          << "<M>/N" << setw(25) <<  "<|M|>/N" << setw(25) \
          << "Xi" << setw(25) << "varM";
  file << "\n";
}

void IsingMonteCarlo::write_exp_vals_to_file(double *expval,ofstream&file, int temp, double varM){
  // write mean energies, magnetization, number of MC cycles,  to file
  file    << setprecision(15) << m_T[temp] << setw(25) << m_MC-m_calibration << setw(25) << m_L2 << setw(25)\
          << expval[0]  << setw(25) <<  expval[1] << setw(25) \
          << expval[2] << setw(25) << varM;
  file << "\n";
}
