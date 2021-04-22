#include "ising_montecarlo.hpp"
#include <random>
#include <chrono>

using namespace std;
using namespace chrono;
void IsingMonteCarlo::wolff_sampling2D(int node_i, int node_j){
  /* Clustering algo (some similarities with breadt-first)

  Input:
  - node_i, node_j: Search keys.  Starting spin of cluster search.
    These input params should be chosen at random
  */
  // -----------------------------------------
   // allocate memory and objects
   //-----------------------------------------

  int sd = chrono::high_resolution_clock::now().time_since_epoch().count() + m_rank;; // Used to obtain seed
  mt19937 gen_accept(sd);                                                       // seeded with sd
  discrete_distribution<> dist_accept({ 1-m_p, m_p});  // creates [0,1] m_p is given in solve.cpp

  int S1, S2, S3, S4;
  // -1 means empty "slot",  assuming no nodes hold negative values (only positive indices)
  int *queue = (int *) malloc((m_L*m_L) * sizeof(*queue));
  int *discovered_i = (int *) malloc((m_L*m_L) * sizeof(*discovered_i));
  int *discovered_j = (int *) malloc((m_L*m_L) * sizeof(*discovered_j));

  // initialize to zeroes and -1, with loop unrolling
  for (int i = 0; i < m_L*m_L; i++){
      discovered_j[i] = 0;
      discovered_i[i] = 0;
      queue[i] = -1;
    }

  int node_id = m_map[node_i]*m_L + m_map[node_j];
  queue[0] = node_id; // index for spin given as search key
  //cluster[node_id] = 1; // store node_id to cluster
  S2d[node_id] *= -1; // flip spin
  discovered_i[0] = node_i;
  discovered_j[0] = node_j;

  bool in_queue = true; // because first queue holds node_id
  // S1,S2 plays this role: int discovered_node; // placeholder to avoid accessing array unnecessary
  int q;               // placeholder for node in queue being worked on
  int cand1, cand2, cand3, cand4;
  // -----------------------------------------
   // Begin algo
   //-----------------------------------------
  while (in_queue == true){
      int i = 0;
      int k = 0;
      //cout << "here" << "\n";
      // enter the adjacent nodes and check if they uphold tau
      while(queue[i] != -1){
        //cout << "q " << queue[i];
        // access the node in queue
        q = queue[i];
        queue[i] = -1; // "delete" the node in queue[i]
        //cout << " " << i << "\n";
        //cout << "check nodes edged to " << q << "\n";
        // check if nodes edged to node (q) in queue
        // is not stored in cluster and upholds criteria to be in cluster

        // four possible neighbours
        cand1 = m_map[discovered_i[i]-1]*m_L + m_map[discovered_j[i]];
        cand2 = m_map[discovered_i[i]+1]*m_L + m_map[discovered_j[i]];
        cand3 = m_map[discovered_i[i]]*m_L + m_map[discovered_j[i]-1];
        cand4 = m_map[discovered_i[i]]*m_L + m_map[discovered_j[i]+1];

        S1 =  S2d[cand1];
        S2 =  S2d[cand2];
        S3 =  S2d[cand3];
        S4 =  S2d[cand4];

        // --> fill in with discovered nodes and flip
        if (S1 != S2d[q]){
            if (dist_accept(gen_accept)){ // draw 0,1 with prob
            queue[k] = cand1; // filling up next qeue with nodes, starting from index 0;
            discovered_i[k] = discovered_i[i]-1;
            discovered_j[k] = discovered_j[i];
            S2d[cand1] *= -1; // flip
            //cout << "add to queue " << discovered_node  << "\n";
            k++;
            }
          }

        if (S2 != S2d[q]){
          if (dist_accept(gen_accept)){
            queue[k] = cand2; // filling up next qeue with nodes, starting from index 0;
            discovered_i[k] = discovered_i[i]+1;
            discovered_j[k] = discovered_j[i];
            S2d[cand2] *= -1; // flip
            //cout << "add to queue " << discovered_node  << "\n";
              k++;
              }
            }

        if (S3 != S2d[q]){
          if (dist_accept(gen_accept)){
            queue[k] = cand3; // filling up next qeue with nodes, starting from index 0;
            discovered_i[k] = discovered_i[i];
            discovered_j[k] = discovered_j[i]-1;
            S2d[cand3] *= -1; // flip
            //cout << "add to queue " << discovered_node  << "\n";
            k++;
          }
        }

        if (S4 != S2d[q]){
            if (dist_accept(gen_accept)){
              queue[k] = cand4; // filling up next qeue with nodes, starting from index 0;
              discovered_i[k] = discovered_i[i];
              discovered_j[k] = discovered_j[i]+1;
              S2d[cand4] *= -1; // flip
              //cout << "add to queue " << discovered_node  << "\n";
              k++;
            }
          }

        //cout << "q " << queue[0] << "\n";
        //cout << "\n";
        i++;
    }
    // checking if there are any nodes in queue
    // after going through queue[0] in the new cycle/depth;
    if (queue[0] == -1){
      in_queue = false; // i.e. end while
    }
  }
  // end algo
}


void IsingMonteCarlo::wolff_sampling1D(int node_id){
  /* Clustering algo (some similarities with breadt-first)

  Input:
  - node_id: Search key.  Starting node of cluster search.
  */

  // -----------------------------------------
   // allocate memory and objects
   //-----------------------------------------

  int sd = chrono::high_resolution_clock::now().time_since_epoch().count() + m_rank;; // Used to obtain seed
  mt19937 gen_accept(sd);                                                       // seeded with sd
  discrete_distribution<> dist_accept({ 1-m_p, m_p});  // creates [0,1] m_p is given in solve.cpp

  // Storing the nodes in cluster dyn. allocate
  // -1 means empty "slot",  assuming no nodes hold negative values
  int *queue = (int *) malloc(m_L* sizeof(*queue));

  // initialize to zeroes and -1, with loop unrolling
  for (int i = 0; i < m_L; i++){
      queue[i] = -1;
    }

  queue[0] = node_id; // put given node as search key
  S2d[node_id] *= -1; // flip spin
  bool in_queue = true; // because first queue holds node_id
  int q;               // placeholder for node in queue being worked on
  int cand1, cand2;
  int S1, S2;
  // -----------------------------------------
   // Begin algo
   //-----------------------------------------
  while (in_queue == true){
      int i = 0;
      int k = 0;
      //cout << "here" << "\n";
      // enter the adjacent nodes and check if they uphold tau
      while(queue[i] != -1){
        //cout << "q " << queue[i];
        // access the node in queue
        q = queue[i];
        queue[i] = -1; // "delete" the node in queue[i]
        //cout << " " << i << "\n";
        //cout << "check nodes edged to " << q << "\n";
        // check if nodes edged to node (q) in queue
        // is not stored in cluster and upholds the SNN criteria

        // two possible neighbours
        cand1 = m_map[q+1]; // to the right
        cand2 = m_map[q-1]; // to the left

        S1 =  S1d[cand1];
        S2 =  S1d[cand2];

        //  access with prob p and flip
        if (S1 != S1d[q]){
            if (dist_accept(gen_accept)){
              queue[k] = cand1; // filling up next qeue with nodes, starting from index 0;
              S2d[cand1] *= -1; // flip
              k++;
            }
          }

        if (S2 != S1d[q]){
              if (dist_accept(gen_accept)){
                queue[k] = cand2; // filling up next qeue with nodes, starting from index 0;
                S1d[cand2] *= -1; // flip
                k++;
              }
            }
        i++;
    }
    // checking if there are any nodes in queue
    // after going through queue[0] in the new cycle/depth;
    if (queue[0] == -1){
      in_queue = false; // i.e. end while
    }
  }
  //end algo
}
