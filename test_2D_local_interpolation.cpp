//This file test the local interpolation sheme for quintic hermite polyomials in two dimension
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>

#include <fstream>
using std::ofstream;
using std::cerr;




#define LSPI_INCLUDE_DEBUG_HELPERS
// #define LSPI_DEBUG_2D
// #define LSPI_DEBUG_1D

#include "LSPI_tools.hpp"
#include "LSPI_2D.hpp"
#include "tictoc.hpp"


// simple data matrix -- row major order
const int gridXsize = 8;
const int gridYsize = 10;
const double gridX[gridXsize] = { 1.0, 2.0, 3.0, 3.7, 4.5, 7.5, 7.75, 8};
const double gridY[gridYsize] = {-1.5,-1.0, 0.0, 1.0, 2.0, 3.0, 3.7 , 4.1, 4.5, 7.5};
const double data[gridYsize*gridXsize] =
      {  -3.0000,   -2.0000,   -1.0000,    1.0000,    3.0000,    2.0000,         0,   -1.0000,
         -2.0000,   -1.0000,    1.0000,    3.0000,    2.0000,         0,   -1.0000,    4.0000,
         -1.0000,    1.0000,    3.0000,    2.0000,         0,   -1.0000,    4.0000,    7.0000,
          1.0000,    3.0000,    2.0000,         0,   -1.0000,    4.0000,    7.0000,    8.0000,
          3.0000,    2.0000,         0,   -1.0000,    4.0000,    7.0000,    8.0000,    9.5000,
          2.0000,         0,   -1.0000,    4.0000,    7.0000,    8.0000,    9.5000,   10.5000,
               0,   -1.0000,    4.0000,    7.0000,    8.0000,    9.5000,   10.5000,   11.0000,
         -1.0000,    4.0000,    7.0000,    8.0000,    9.5000,   10.5000,   11.0000,   -3.0000,
          4.0000,    7.0000,    8.0000,    9.5000,   10.5000,   11.0000,   -3.0000,    2.0000,
          7.0000,    8.0000,    9.5000,   10.5000,   11.0000,   -3.0000,         0,         0
      };
// const double data[gridYsize*gridXsize] = {
//     1.1,  1.2,  1.3,  1.4,  1.5,  1.6,  1.7,  1.8,
//     2.1,  2.2,  2.3,  2.4,  2.5,  2.6,  2.7,  2.8,
//     3.1,  3.2,  3.3,  3.4,  3.5,  3.6,  3.7,  3.8,
//     4.1,  4.2,  4.3,  4.4,  4.5,  4.6,  4.7,  4.8,
//     5.1,  5.2,  5.3,  5.4,  5.5,  5.6,  5.7,  5.8,
//     6.1,  6.2,  6.3,  6.4,  6.5,  6.6,  6.7,  6.8,
//     7.1,  7.2,  7.3,  7.4,  7.5,  7.6,  7.7,  7.8,
//     8.1,  8.2,  8.3,  8.4,  8.5,  8.6,  8.7,  8.8,
//     9.1,  9.2,  9.3,  9.4,  9.5,  9.6,  9.7,  9.8,
//    10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8};

using namespace LSPI;





/**
 * \brief Returns the interpolation result in a freshly new()ed array. Remember to delete[] it!
 */
double* runtest(int nrep, const double* qx, const double* qy, int qx_size, int qy_size, int derivatives) {
    displayText("Running test for derivatives = %d ... ", derivatives);
    // caching: create work spaces for interpolation and index pointers for remembering the last evaluation patch
    const int ws = 8;   // workspace size
    double ipparam_xz[ws], ipparam_xzh[ws], ipparam_yz[ws], ipparam_yzh[ws];
    int idx_old_x = -1;
    int idx_old_y = -1;
    // allocate result array
    double* values = new double[qy_size * qx_size]; // NEWed storage for interpolated values
    // perform repeated interpolation
    int ticID = tictoc::tic();
    for (int n=0; n<nrep; ++n) {
        for (int i=0; i<qy_size; ++i){
            for (int j=0; j<qx_size; ++j){ 
                // printf("i = %2d  j = %2d  ", i, j)
                double value = LSPI::compute_interpolation_2D( 
                    qx[j], qy[i], gridX, gridY, data, gridXsize, gridYsize, derivatives, 0.75, 1.0,
                    ipparam_xz, ipparam_xzh, idx_old_x, ipparam_yz, ipparam_yzh, idx_old_y);
                // printf("val = %10.3f\n", value);
                LSPI_MATRIX_ELEMENT(values, i, j, qx_size) = value;
            }
        }
    }
    double duration = tictoc::toc(ticID);
    displayText(" --> Took %8g ms.\n", duration);
    return values;
}


bool arrays_identical(const double* ar1, const double* ar2, size_t length) {
    for (size_t i=0; i < length; ++i)
        if (ar1[i] != ar2[i]) return false;
    return true;
}



int main(int argc, char** argv)
{
    // number of repetitions
    int nrep      = 1;
    int numpoints = 10;
    int output    = 0;

    // Display short help upon "--help"
    if (argc >= 2 && std::string(argv[1]) == "--help") {
        std::cout << "Usage:  " << argv[0] << " silent n nrep" << std::endl;
        std::cout << "      n: number of evaluation points between grid points [default if not specified: " << numpoints << "]" << std::endl;
        std::cout << "   nrep: number of repetitions                           [default if not specified: " << nrep      << "]" << std::endl;
        std::cout << " silent: 0: display only runtime --- 1: display interpolated values (Matlab format);"                     << std::endl;
        return 0;
    }

    if (argc >= 2 && argv[1] != nullptr) numpoints = std::stoi(argv[1]);
    if (argc >= 3 && argv[2] != nullptr)      nrep = std::stol(argv[2]);
    if (argc >= 4 && argv[3] != nullptr)    output = std::stol(argv[3]);

    // display config
    std::cout << "numpoints = " << numpoints << std::endl;
    std::cout << "     nrep = " << nrep      << std::endl;
    std::cout << "   output = " << output    << std::endl;

    // setup query points between grid points
    double  a = 0;
    double  b = 1;
    double* linp = new_linspaced_darray(a,b,numpoints);
    
    // dimensions of query points
    const int qxsize = (gridXsize-1)*(numpoints-1)+1;
    const int qysize = (gridYsize-1)*(numpoints-1)+1;

    // create query points
    double qx[static_cast<size_t>(qxsize)];
    double qy[static_cast<size_t>(qysize)];
    qx[0] = gridX[0];
    qy[0] = gridY[0];
    int idx = 1;
    double griddiff;
    // query point x values
    for(int j=0; j<(gridXsize-1); ++j){
        griddiff = gridX[j+1] - gridX[j];
        for(int i=1; i<=(numpoints-1); ++i){
            qx[idx++] = gridX[j] + linp[i] * griddiff;
        }
    }
    // query point y values
    idx = 1;
    for(int j=0; j<(gridYsize-1); ++j){
        griddiff = gridY[j+1] - gridY[j];
        for(int i=1; i<=(numpoints-1); ++i){
            qy[idx++] = gridY[j] + linp[i] * griddiff;
        }
    }

    // perform interpolation
    double* values0 = runtest(nrep, qx, qy, qxsize, qysize, 0);
    double* values1 = runtest(nrep, qx, qy, qxsize, qysize, 1);
    double* values2 = runtest(nrep, qx, qy, qxsize, qysize, 2);

    // display 
    if (output > 0) {
        // display data
        displayVec(gridX, gridXsize, "x");
        displayVec(gridY, gridYsize, "y");
        displayMat(data, gridYsize, gridXsize, "D");
        // display results for matlab
        displayVec(qx, qxsize, "qx");
        displayVec(qy, qysize, "qy");
        displayMat(values0, qysize, qxsize, "D0");
        displayMat(values1, qysize, qxsize, "D1");
        displayMat(values2, qysize, qxsize, "D2");
    }

    // single point test
    double qqx[1]{ 3.7 };
    double qqy[1]{ 2.0 };
    displayText("Testing single point evaluation\n");
    double* sval = runtest(1, qqx, qqy, 1, 1, 0);
    std::cout << "xquery = " << *qqx  << std::endl;
    std::cout << "yquery = " << *qqy  << std::endl;
    std::cout << "xyval  = " << *sval << std::endl;
    delete[] sval;

    // test on data: same values returned?
    displayText("Checking interpolation at grid points ...\n");
    double* values = runtest(1, gridX, gridY, gridXsize, gridYsize, 0);
    if (arrays_identical(values, data, gridXsize*gridYsize)) {
        displayText("OK (identical)\n");
    } else {
        displayText("FAILURE -- interpolation property violated\n");
    }

    if (output > 0) {
        displayVec(linp , numpoints, "linp");
        displayVec(gridX, gridXsize, "gridX");
        displayVec(gridY, gridYsize, "gridY");
        displayVec(data , gridYsize*gridXsize, "data1D");
        displayMat(data , gridYsize,gridXsize, "data2D");
        displayMat(values, gridYsize, gridXsize, "interpolation2D");
    }

    // cleanup
    delete[] linp;
    delete[] values0;
    delete[] values1;
    delete[] values2;

    return 0;
}
