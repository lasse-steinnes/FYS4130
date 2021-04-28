# FYS4130
Git for Assignment 2 in the UiO course FYS4130 (Statistical Mechanics)

### Main overview
* The programs in this repository aim at finding expectation values for a 1D and 2D Ising Model, studying phase transition in a ferromagnetic system.

* The main challenge was to apply Monte Carlo simulations and the Wolff sampling algorithm to obtain expectation values and study phase transition for a 1D and 2D system. To do so a random number generator was used (see link below). In addition the algorithm was parallelized (over temperatures) using OpenMP. Each thread was regarded as a separate experiment, thus each thread has its own unique seed. The spin matrix is initialized once for each thread. For temperatures within that thread, the first Monte Carlo (MC) cycle of the next temperature uses the last spin configuration of the previous temperature.

* Textfiles and figures can be found in the folder Results.

* The calculations are performed with the energy coefficent J, and boltzmann factor, k, set to 1. With this scaling the critical temperature is approximately 2.269 in the termodynamical limit L -> infinity (infinite number of spins in 2D system). Termodynamical properties are scaled by number of spins (1D: L/ 2D: L^2), since the expectation values are an extensive parameters, meaning that the amplitude of magnetic moment depends on the size of the 1D/2D spin system.

### Code: Link and description of programmes
- [initalizer.cpp](https://github.com/lasse-steinnes/FYS4130/blob/main/assignment2/initializer.cpp) : Initializes parameters and spin vectors.

- [main.cpp](https://github.com/lasse-steinnes/FYS4130/blob/main/assignment2/main.cpp) : Runs the other programmes and provide user options through terminal.

- [magnetic_moment.cpp](https://github.com/lasse-steinnes/FYS4130/blob/main/assignment2/magnetic_moment.cpp) : Calculates the total magnetic moment of either a 1D chain or 2D grid of spins.

 - [makefile](https://github.com/lasse-steinnes/FYS4130/blob/main/assignment2/makefile) : Compiles and executes cpp files.

-  [ising_montecarlo.hpp](https://github.com/lasse-steinnes/FYS4130/blob/main/assignment2/ising_montecarlo.hpp) : Headerfile for the IsingMonteCarlo class.

- [solve.cpp](https://github.com/lasse-steinnes/FYS4130/blob/main/assignment2/solve.cpp) : Performs the MC simulation with Wolff sampling.

- [wolf_sample.cpp](https://github.com/lasse-steinnes/FYS4130/blob/main/assignment2/wolf_sample.cpp) : Provdes Wolff sampling for the 1D and 2D case.

- [random_draw.cpp](https://github.com/lasse-steinnes/FYS4130/blob/main/assignment2/random_draw.cpp) : Helps initialize the system.

- [write.cpp](https://github.com/lasse-steinnes/FYS4130/blob/main/assignment2/write.cpp) : write to file code.

- [plot_corr.py](https://github.com/lasse-steinnes/FYS4130/blob/main/assignment2/plot_corr.py) plots the correlation function (1D case).

- [plot_gamma.py](https://github.com/lasse-steinnes/FYS4130/blob/main/assignment2/plot_gamma.py) plots the gamma function to obtain the critical temperature for the 2D case.

- [plotm_and_m2.py](https://github.com/lasse-steinnes/FYS4130/blob/main/assignment2/plotm_and_m2.py) plots <m> and <m^2>.



The cpp-files can be compiled in terminal with "make all", and can be run by ./main.out.

### Links and packages
- The Mersenne Twister (pseudo)random number generator was used in generating uniform distribution to draw indices and acceptance criteria. Documentation on the class mt19937_64 can be found [here.](https://www.cplusplus.com/reference/random/mt19937_64/)

- Documentation for Matplotlib from python from [here](https://matplotlib.org/)

- Documentation on parallelization with OpenMP can be found [here](https://www.openmp.org/wp-content/uploads/OpenMP-4.5-1115-CPP-web.pdf) or for more versions [here](https://www.openmp.org/resources/refguides/)
