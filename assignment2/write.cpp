// Programmes to write expectation values to file
#include "ising_montecarlo.hpp"


void IsingMonteCarlo::open_exp_vals_to_file(ofstream&file){ // write expectation values after all cycles
  string filename("./Results/Data/expvalues" + to_string(m_MC) + \
                  "-" + to_string(m_L) + "by" + to_string(m_L) + "rank" + to_string(m_rank) + ".txt");
  file.open(filename);
  file    << "T" << setw(25) << "MC_cycles-ac" << setw(25) << "N_spins" << setw(25)\
          << "<M>/N" << setw(25) <<  "<M^2>/N" << setw(25) \
          << "<M^4>" << setw(25) << "Gamma";
  file << "\n";
}

void IsingMonteCarlo::write_exp_vals_to_file(double *expval,ofstream&file, double temperature, double gamma){
  // write mean energies, magnetization, number of MC cycles,  to file
  file    << setprecision(15) << temperature << setw(25) << m_MC << setw(25) << m_L2 << setw(25)\
          << expval[0]  << setw(25) <<  expval[1] << setw(25) \
          << expval[2] << setw(25) << gamma;
  file << "\n";
}

void IsingMonteCarlo::write_corr(double *arr, int temp){
  // open spin to file if true
  ofstream spinfile;
  string filename("./Results/Data/Corr1D" + to_string(m_MC) + \
                  "-" + to_string(m_L) + "by" + to_string(m_L) + "rank" + to_string(m_rank) + ".txt");
  spinfile.open(filename);

  // write spin file
  spinfile << setprecision(4) <<  "T:" << " " << m_T[temp] << "\n";
  for (int i = 0; i < m_L; i++){
      spinfile << setprecision(15) <<  arr[i]; // get spin matrix
      spinfile << "\n";
      }
  // close spin file
  spinfile.close(); //
  }
