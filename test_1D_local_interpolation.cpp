//This file test the local interpolation sheme for quintic hermite polyomials
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>
#include <random>
#include <algorithm>


#define LSPI_INCLUDE_DEBUG_HELPERS
// #define LSPI_DEBUG_1D

#include "LSPI_tools.hpp"
#include "LSPI_1D.hpp"
#include "tictoc.hpp"

//////////////////////////////////////////////////////////////////////


using namespace LSPI;


/**
 * This function runs the interpolation and times it
 */
double runtest(int nrep, int derivatives, const double* pvecIN, int plen, const double* grid, int gridlen, const double* data, double alpha, double* values, unsigned int randseed) {
   std::cout<<"Doing "<<nrep<<" repetitions with derivatives = "<<derivatives<<" -- ";
   int ticID;                                  // carrier for time measurement
   double duration_random = 0;                 // runtime for random number generations
   double duration_interp = 0;                 // runtime for interpolation
   int idx_old = -1;                           // for caching: last found interval
   double interp_param[8]{};                   //              and respective coefficients
   double* pvec = new double[plen];            // local copy of input data
   std::copy(pvecIN, pvecIN + plen, pvec);
   // randomize evaluation positions - if requested
   if (randseed > 0){
      ticID = tictoc::tic();
      std::mt19937 randgen(randseed);
      std::shuffle(pvec, pvec + plen, randgen);
      duration_random += tictoc::toc(ticID);
   }
   // interpolation
   ticID = tictoc::tic();
   for (int r=0; r<nrep; ++r) {
      for(int i=0; i<plen; ++i){
         values[i] = compute_interpolation_1D(pvec[i], grid, data, gridlen, derivatives, alpha, interp_param, idx_old);
      }
   }
   duration_interp = tictoc::toc(ticID);
   std::cout << std::endl;
   std::cout << "Total runtime: " << duration_interp + duration_random << " ms" << std::endl;
   std::cout << "Interpolation: " << duration_interp << " ms" << std::endl;
   std::cout << "Randomization: " << duration_random << " ms" << std::endl;
   delete[] pvec;
   return duration_interp + duration_random;
}



int main(int argc, char** argv)
{
   // Defaults
   int numpoints = 10;
   int nrep = 1000000;
   int output = 0;
   unsigned int randseed = 0;
   
   // Display short help upon "--help"
   if (argc >= 2 && std::string(argv[1]) == "--help") {
      std::cout << "Usage:  " << argv[0] << " n nrep output randseed" << std::endl;
      std::cout << "          n: number of evaluation points between grid points  [default if not specified: " << numpoints << "]" << std::endl;
      std::cout << "       nrep: number of repetitions for time measurement       [default if not specified: " <<      nrep << "]" << std::endl;
      std::cout << "     output: 0: display  runtime  -- 1: display also values   [default if not specified: " <<    output << "]" << std::endl;
      std::cout << "   randseed: random seed for shuffling evaluation points      [default value 0 disables randomization]"        << std::endl;
      return 0;
   }

   // Hard coded example data  //  0                                      1                                                   2
   const int gridlen = 28;     //  1    2    3   4   5   6   7    8   9   0    1   2     3    4     5    6    7     8    9    0    1    2    3    4    5    6    7    8
   // this example uses non-equidistant grid -- so no guaranteed runtime
   const double grid[gridlen] ={-1.5,-1.0, 0.0,1.0,2.0,3.0,3.7, 4.1,4.5,7.5,7.75,8.0, 8.25, 8.5, 8.75, 9.0,10.0,10.25,10.5,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,21.0};
   const double data[gridlen] ={-3.0,-2.0,-1.0,1.0,3.0,2.0,0.0,-1.0,4.0,7.0,8.0 ,9.5,10.0 ,11.0,-3.0 ,15.0, 1.0,10.0 , 9.0, 8.0, 7.0, 6.0, 5.0, 3.5, 3.0, 2.0, 1.0, 0.0};
   std::cout<<"grid with "; LSPI::print_vector(grid,gridlen);
   std::cout<<"data with "; LSPI::print_vector(data,gridlen);

   // Retrieve command line args - no error handling
   if (argc >= 2 && argv[1] != nullptr) numpoints = std::stoi(argv[1]);
   if (argc >= 3 && argv[2] != nullptr)      nrep = std::stol(argv[2]);
   if (argc >= 4 && argv[3] != nullptr)    output = std::stol(argv[3]);
   if (argc >= 5 && argv[4] != nullptr)  randseed = std::stoul(argv[4]);  // ignore overflow

   // Display config
   std::cout << "numpoints = " << numpoints << std::endl;
   std::cout << "     nrep = " << nrep      << std::endl;
   std::cout << "   output = " << output    << std::endl;
   std::cout << " randseed = " << randseed  << std::endl;


   // Create helper for spacing between grid points for evaluation
   const double a = 0;
   const double b = 1;
   const double* lin01 = new_linspaced_darray(a,b,numpoints);

   // Create evaluation points
   const int plen = (gridlen-1)*(numpoints-1)+1;
   double pvec[static_cast<size_t>(plen)];  // target vector of points to interpolate at
   pvec[0] = grid[0];  // initial point is first grid point
   int idx = 1;      
   for(int j=0;j<(gridlen-1);j++){
      double difference = grid[j+1]-grid[j];
      for(int i=1;i<=(numpoints-1);i++){
         // std::cout<< "idx: "<<idx<< ", i: "<<i << ", j: "<<j<< ", i+j:"<< i+j << ", and lin01[i]" << lin01[i]<<std::endl;
         pvec[idx] = grid[j] + lin01[i] * difference;
         idx++;
      }
   }
   delete[] lin01;

   // initialize output variables
   double valuesD0[static_cast<size_t>(plen)], durationD0;
   double valuesD1[static_cast<size_t>(plen)], durationD1;
   double valuesD2[static_cast<size_t>(plen)], durationD2;

   durationD0 = runtest(nrep, 0, pvec, plen, grid, gridlen, data, 0.3, valuesD0, randseed);
   durationD1 = runtest(nrep, 1, pvec, plen, grid, gridlen, data, 1.0, valuesD1, randseed);
   durationD2 = runtest(nrep, 2, pvec, plen, grid, gridlen, data, 1.0, valuesD2, randseed);

   std::cout << "%% D0: duration = " << durationD0 << std::endl;
   std::cout << "%% D1: duration = " << durationD1 << std::endl;
   std::cout << "%% D2: duration = " << durationD2 << std::endl;

   if (output > 0) {
      std::cout << "%% DO values: ";  LSPI::print_vector(valuesD0,plen);  std::cout << std::endl;
      std::cout << "%% D1 values: ";  LSPI::print_vector(valuesD1,plen);  std::cout << std::endl;
      std::cout << "%% D2 values: ";  LSPI::print_vector(valuesD2,plen);  std::cout << std::endl;
   }

} // end of main



