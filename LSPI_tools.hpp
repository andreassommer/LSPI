#ifndef LSPI_TOOLS_HEADER
#define LSPI_TOOLS_HEADER

/**
 *  \file Collection of tools for local shape preserving interpolation 
 *  \author Andreas Sommer <code@andreas-sommer.eu>
 *  \author Julian Niederer <julian.niederer@iwr.uni-heidelberg.de>
 *  \date   Jul2026
 */

// includes required for LSPI_1D and LSPI_2D
#include <cstring>
#include <cmath>

/** 
 * DEBUG MACROS -- verbosity
 */

#ifdef LSPI_DEBUG_1D
   #define LSPI_VERBOSE_1D(x) x
#else
   #define LSPI_VERBOSE_1D(x)
#endif

#ifdef LSPI_DEBUG_2D
   #define LSPI_VERBOSE_2D(x) x
#else
   #define LSPI_VERBOSE_2D(x)
#endif

#if defined(LSPI_DEBUG_1D) || defined(LSPI_DEBUG_2D)
   #define LSPI_INCLUDE_DEBUG_HELPERS
#endif

#if defined(LSPI_INCLUDE_DEBUG_HELPERS)
   #include <iostream>
   #include <cstdlib>
   #include <cstdarg>
#endif


/**
 * \brief Matrix accessor: row-major order
 */
#define LSPI_MATRIX_ELEMENT(matrix, row, col, ncols)  (matrix[(row)*(ncols)+(col)])



namespace LSPI {


const char* const DEFAULT_NUMBER_FORMAT = " %12.6g ";   // keep spaces! otherwise: -4.111-3.222 may happen!
const char*       NUMBERFORMAT = DEFAULT_NUMBER_FORMAT; // NUMBERFORMAT can point to different formats


/**
 * \brief Enum to specify X or Y direction for 2D case
 */
enum Direction {
    X, Y
};


/**
 * \brief Copies parts of source into array 2
 * \param[in] source source array
 * \param[in] idx1   start index for source
 * \param[in] count  number of elements to be copied
 * \param[in] target target array
 * \param[in] idx2   start index for target
 */
template<typename T>
inline
void copy_entries(const T* source, int idx1, int count, T* target, int idx2)
{
   for(int i=0; i<count; ++i){
      target[idx2+i] = source[idx1+i];
   }
}



/**
 * \brief Calculate 5th order smoothstep function.
 * \param[in]  u: position of evaluation, scaled to [0,1]
 * \return Value of polynomial: -6u^5 + 15u^4 - 10u^3 + 1
 */
template <typename T> 
inline
T smoothstep5(const T &u)
{
    return u * u * u * ( u * (-6.0 * u + 15.0) - 10.0 ) + 1.0;
}



/**
 * \brief Get the signum as one or minus-one in the respective type
 */
template <typename T>
inline
T get_signum(T y)
{
   if (y>=0)
      return (T) 1.0;
   else
      return (T) -1.0;
}



/**
 * \brief Computes elementwise difference of inputArray and store result in diffArray
 * \param[in]  inputArray array to compute element-wise differences from
 * \param[in]  size       size of inputArray
 * \param[out] diffArray  array with element-wise differences from inputArray
 */
template <typename T>
inline
void diff(const T* inputArray, int length, T* diffArray)
{
   for (int i = 0; i < length; ++i) {
      diffArray[i] = inputArray[i+1] - inputArray[i];
   }
}



/**
 * \brief Computes quotient of two arrays at specific index  [cpq = compute pointwise quotient]
 * \param[in]  inputArray array to compute element-wise differences from
 * \param[in]  size       size of inputArray
 * \param[out] diffArray  array with element-wise differences from inputArray
 */
template <typename T>
inline
T cpq(const T* array1, const T* array2, int idx)
{
   return array1[idx] / array2[idx];
}



/**
 * \brief Computes element-wise quotients of two arrays and stores result in array [cpq = compute pointwise quotient]
 * \param[in]  numerator      array to compute element-wise differences from
 * \param[in]  denominator    array to compute element-wise differences from
 * \param[in]  size           size of input arrays
 * \param[out] quotientArray  array with element-wise quotients 
 */
template <typename T>
inline
void cpq_array(const T* numerator, const T* denominator, int length, T* quotientArray)
{
   for(int i = 0; i < length; ++i) {
      quotientArray[i] = numerator[i] / denominator[i];
   }
}



/**
 * \brief Creates new linspaced double vector (for testers)
 */
inline
double* new_linspaced_darray(double start, double end, int count)
{
   double* linspacedvector = new double[count];
   double delta = (end - start) / (count - 1);
   for(int i=0; i < count-1; ++i)
   {
      linspacedvector[i] = (start + delta * i);
   }
   linspacedvector[count-1] = end;
   return linspacedvector;
}




/**
 * \brief Find the interval in grid in where point p resides
 * \param[in] p         query point
 * \param[in] grid      grid array - assumed to be strictly monotonically increasing
 * \param[in] gridsize  size of grid array
 * \return Interval i where p resides: grid[i] <= p < grid[i+1]
 * \note See find_index_from() for a fast version if an old index is available
 */
template <typename T>
int find_index(const T& p, const T* grid, int gridsize)
{
   // SPEEDUP POSSIBLE: start search from idxguess forward or backward --> see find_index_from()
   int idx = 0;
   double step = grid[1]-grid[0];  // lucky guess if equidistant
   int idxguess = static_cast<int>(p / step);
   if( idxguess >= gridsize){idxguess = 0;}
   if (idxguess < 0) { idxguess = 0; }
   if (grid[idxguess] > p){idxguess = 0;} 
   // find index of interval, where x lies in 
   for (idx = idxguess; idx < gridsize; ++idx) {
      if (p < grid[idx]) break;
   }
   // ensure 
   if(idx >= gridsize){
      idx = idx -2;
   }else{
      idx = idx - 1;
   }
   if(idx < 0){
      idx = 0;
   }
   return idx;
}



/**
 * \brief Find the interval in grid in where point p resides - starting from given index
 * \param[in] p         query point
 * \param[in] grid      grid array - assumed to be strictly monotonically increasing
 * \param[in] gridsize  size of grid array
 * \param[in] startidx  index at where to start the search
 * \return Interval i where p resides: grid[i] <= p < grid[i+1]
 */
template <typename T>
int find_index_from(const T& p, const T* grid, int gridsize, int startidx)
{
   if (startidx < 0 || startidx >= gridsize)
      return find_index(p, grid, gridsize);          // fallback on bad input
   int idx = 0;
   if (p <= grid[startidx]){                         // p <= grid[startidx]
      for (idx = startidx; idx >= 0; --idx)          // 
         if (grid[idx] <= p) break;
   }else{                                            // p  > grid[startidx]
      for (idx = startidx+1; idx < gridsize; ++idx)  // start at next index and search for larger grid point
         if (grid[idx] > p) { idx--; break; }        // have to go back one because we check for > p
   }
   // safety code
   if (idx < 0) idx = 0;
   if (idx >= gridsize) idx = gridsize-2;            // have to go back 2: last interval starts at that node
   return idx;
}




// ******************************************************************************************
// ********* DEBUG HELPERS ******************************************************************
#ifdef LSPI_INCLUDE_DEBUG_HELPERS

template<typename T>
void displayVec(const T* vec, int vecsize, const char* name = "vec") {
   const char* fmt = LSPI::NUMBERFORMAT;
   if (name)  printf("%s = ", name);
   printf("[ ");
   for (int i = 0; i < vecsize; ++i) {
      printf(fmt, vec[i]);
   }
   printf(" ];\n");
}

template<typename T>
void displayMat(const T* matrix, int nrows, int ncols, const char* name = "vec") {
   const char* fmt = LSPI::NUMBERFORMAT;
   if (name)  printf("%s = ", name);
   printf("[ \n");
   for (int row = 0; row < nrows; ++row) {
      for (int col = 0; col < ncols; ++col) {
         printf(fmt, LSPI_MATRIX_ELEMENT(matrix, row, col, ncols));
      }
      if (row < nrows-1)  printf(" ;\n");
   }
   printf(" ];\n");
}

template<typename T>
void displayVar(const T& var, const char* name = "var", const char* separator = "\n") {
   const char* fmt = LSPI::NUMBERFORMAT;
   if (name)
      printf("%s = ", name);
   printf(fmt, (double) var);
   if (separator)
      printf(separator);
}

int displayFmt(const char* fmt, ...){
    va_list args;
    va_start(args, fmt);
    int rc = std::vprintf(fmt, args);
    va_end(args);
    return rc;
}

int displayText(const char* fmt, ...){
    printf("%% ");
    va_list args;
    va_start(args, fmt);
    int rc = std::vprintf(fmt, args);
    va_end(args);
    return rc < 0 ? rc : rc + 2; // +2 for "% "
}

/**
 * \brief Debug helper: print a vector to stdout
 */
template<typename T>
inline
void print_vector(const T* vec, int length)
{
std::cout << "given length: " << length << std::endl<<"[";
for (int i=0; i<length;++i)
   std::cout << vec[i] << ", ";
std::cout <<"]"<< std::endl;
}

template<typename T>
void print_vector_2D(const T* data,const int &length_y, const int &length_x)
{
  //we assume one dim vector which is interpreted as a two dim array.  
  std::cout << "Matrix with given length: " << length_x*length_y << std::endl<< "[";
  for (int i=0; i<length_y;++i){
    for(int j=0; j<length_x;++j){
        if (j != 0)  std::cout << ", ";
        std::cout << LSPI_MATRIX_ELEMENT(data, i, j, length_x);
    }
    std::cout<<";"<<std::endl;
  }
  std::cout<<"]"<<std::endl;
}

#endif // end of LSPI_INCLUDE_DEBUG_HELPERS
// ******************************************************************************************
// ******************************************************************************************






} // end of namespace

#endif // include guard