#include "ising_montecarlo.hpp"
#include <random>
#include <chrono>

using namespace std;
using namespace chrono;

/* 2D code */

void IsingMonteCarlo::growCluster2D(int i, int j, int clusterSpin) {

    // mark the spin as belonging to the cluster and flip it
    cluster[i*m_L+j] = true;
    S2d[i*m_L +j] = -1*S2d[i*m_L + j];
    m_MagneticMoment += 2*S2d[i*m_L + j]; // change magnetic moment

    // find the indices of the 4 neighbors
    // assuming periodic boundary conditions
    int iPrev = i == 0    ? m_L-1 : i-1;
    int iNext = i == m_L-1 ? 0    : i+1;
    int jPrev = j == 0    ? m_L-1 : j-1;
    int jNext = j == m_L-1 ? 0    : j+1;

    // if the neighbor spin does not belong to the
    // cluster, then try to add it to the cluster
    if (!cluster[iPrev*m_L + j])
        tryAdd2D(iPrev, j, clusterSpin);
    if (!cluster[iNext*m_L + j])
        tryAdd2D(iNext, j, clusterSpin);
    if (!cluster[i*m_L + jPrev])
        tryAdd2D(i, jPrev, clusterSpin);
    if (!cluster[i*m_L + jNext])
        tryAdd2D(i, jNext, clusterSpin);
}

void IsingMonteCarlo::tryAdd2D(int i, int j, int clusterSpin){
    if (S2d[i*m_L + j] == clusterSpin){
        if (m_distribution(m_gen) < m_p){
            growCluster2D(i, j, clusterSpin);
          }
        }
}


/* 1D code */

void IsingMonteCarlo::growCluster1D(int i, int clusterSpin) {
    // mark the spin as belonging to the cluster and flip it
    cluster1D[i] = true;
    S1d[i] = -1*S1d[i];
    m_MagneticMoment += 2*S1d[i]; // change magnetic moment

    // find the indices of the 4 neighbors
    // assuming periodic boundary conditions
    int iPrev = i == 0    ? m_L-1 : i-1;
    int iNext = i == m_L-1 ? 0    : i+1;

    // if the neighbor spin does not belong to the
    // cluster, then try to add it to the cluster
    if (!cluster1D[iPrev])
        tryAdd1D(iPrev, clusterSpin);
    if (!cluster1D[iNext])
        tryAdd1D(iNext, clusterSpin);
}

void IsingMonteCarlo::tryAdd1D(int i, int clusterSpin){
    if (S1d[i] == clusterSpin){
        if (m_distribution(m_gen) < m_p){
            growCluster1D(i, clusterSpin);
          }
        }
}
