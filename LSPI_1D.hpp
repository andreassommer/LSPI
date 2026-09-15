#ifndef LSPI_1D_HEADER
#define LSPI_1D_HEADER

/**
 *  \file   Header for 1D local shape preserving interpolation 
 *  \author Julian Niederer <julian.niederer@iwr.uni-heidelberg.de>
 *  \author Andreas Sommer <code@andreas-sommer.eu>
 */


#include "LSPI_tools.hpp"


namespace LSPI {



/**
 * \brief Compute weighted harmonic mean (whm) or decide to use 0.
 * \input[in]
 */
template <typename T>
T compute_weighted_harmonic_mean(const T& delta_im, const T& delta_i,
                                 const T& h_im    , const T& h_i,
                                 const T& hdm_im  , const T& hdm_i,
                                 const T& hdp_im  , const T& hdp_i,
                                 const T& gamma   , const T& left)
{
   LSPI_VERBOSE_1D(
      displayText("info for whm:  hdm_i = %g  hdp_im 0 %g ", hdm_i, hdp_im);
      bool c1 = (delta_im*delta_i>0);   // c1: check monotonicity and ensure not to divide by zero
      bool c2 = (left>=hdm_im+hdp_im);  // c2: checks on the grid to get a possible solution
      bool c3 = (left>=hdm_i +hdp_i);   // c3: checks on the grid to get a possible solution
      bool c4 = (left>=2*hdm_i);        // c4: check sufficient condition for monoticity .. ,    // 3 - 2 * h_i/h_i-1 >=0
      bool c5 = (left>=2*hdp_im);       // c5: .. preservation when "rule" in STEP 3 is applied. // 3 - 2 * h_i-1/h_i >=0
      displayText("Conditions: %d %d %d %d %d\n", c1, c2, c3, c4, c5);
   )
   if (delta_im*delta_i <= 0)  return (T) 0.0;   // condition c1
   if (left < hdm_im+hdp_im)   return (T) 0.0;   // condition c2
   if (left < hdm_i +hdp_i)    return (T) 0.0;   // condition c3
   if (left < 2*hdm_i)         return (T) 0.0;   // condition c4
   if (left < 2*hdp_im)        return (T) 0.0;   // condition c5
   // here, all checks are passed and we compute the weighted harmonic mean
   return (T) gamma * (h_im  + h_i) / (h_i/delta_im  + h_im /delta_i);
}



/**
 * \brief Evaluate weighted hermite basis functions -- with NO DERIVATIVE weights
 * \param[in] u   point of evaluation in the respective interval (p-x_i)/h_i
 * \param[in] interp_param = [g,h_i,z_i,z_i+1,p_i,p_ip] 
 */
template <typename T>
inline
T eval_weighted_hermite_basis_functions_DERIV0(const T& u, const T* interp_param){
   T mu  = 1-u;
   T u2  = u*u;
   T umu = u*mu;
   T row = umu *(1+u-9*u2+6*u2*u);
   return (T) ( interp_param[2]*(mu + row) 
              + interp_param[3]*(u  - row));
}



/**
 * \brief Evaluate weighted hermite basis function - only FIRST derivative weights p_i, p_ip
 * \param[in] u: point of evaluation in the respective interval (p-x_i)/h_i
 * \param[in] interp_param = [g,h_i,z[2],z[3],p_i,p_ip] 
 * \return value of hermite basis function
 */
template <typename T>
inline
T eval_weighted_hermite_basis_functions_DERIV1(const T& u, const T* interp_param){
   T mu  = 1-u;
   T mu2 = mu*mu;
   T u2  = u*u;
   T umu = u*mu;
   T row = umu*(1+u-9*u2+6*u2*u);
   return ( interp_param[2]*(mu + row)
          + interp_param[3]*(u  - row)
          + interp_param[4]*(umu*interp_param[1]*mu2*(3*u+1))
          + interp_param[5]*(umu*interp_param[1]* u2*(3*u-4)) );
}



/**
 * \brief Evaluate weighted hermite basis function - with FIRST and SECOND derivative weights p_i, p_ip, Pi, P_ip
 * \param[in] u: point of evaluation in the respective interval (p-x_i)/h_i
 * \param[in] interp_param = [g,h_i,z[2],z[3],p_i,p_ip,P_i,P_ip]
 * \return value of hermite basis function
 */
template <typename T>
T eval_weighted_hermite_basis_functions_DERIV2(const T& u, const T* interp_param){
   T mu = 1-u;
   T mu2 = mu*mu;
   T u2 = u*u;
   T umu = u*mu;
   T umuh = umu*interp_param[1];
   T row = umu*(1+u-9*u2+6*u2*u);
   return (   interp_param[2]*( mu + row)
            + interp_param[3]*(  u - row)
            + interp_param[4]*(umu *interp_param[1]*mu2*(3*u+1))
            + interp_param[5]*(umu *interp_param[1]*u2 *(3*u-4))
            + interp_param[6]*(umuh*interp_param[1]*u  *mu2*0.5)
            + interp_param[7]*(umuh*interp_param[1]*u2 *mu *0.5));
}



/**
 * \brief Calculate coefficients for the weighted harmonic mean.
 * \param[in] p       the harmonic mean at position i
 * \param[in] h_im    h in interval i-1 
 * \param[in] h_i     h in interval i
 * \param[in] sign_im in interval i-1
 * \param[in] sign_i  in interval i
 * \param[in] gamma   coefficient \in (0,1]
 * \return Returns p if acceptable for second derivative, 0.0 if not.
 */
template <typename T>
T calc_P(const T& p, const T& h_im, const T& h_i, const T& sign_im, const T& sign_i, const T& gamma){
   if ( p == 0.0 ){
      LSPI_VERBOSE_1D( displayText("Conditions for P: NONE, already zero\n"); )
      return (T) 0.0;
   }else{
      T cond_p = 4*p;
      bool c1 = (20.0-gamma*(8+h_i)*(1+h_i/h_im)>=0);
      bool c3 = (        sign_im*(cond_p-h_im*p)>=0);
      bool c2 = (        sign_i *(cond_p+h_i *p)>=0);
      LSPI_VERBOSE_1D( displayText("Conditions for P: %d %d %d\n", c1, c2, c3); )
      if( c1 && c2 && c3){
         return (T) p;
      }else{
         return (T) 0.0;
      }
   }
}



/**
 * \brief Fills a window of `count` points around `start`, extrapolating linearly beyond the array boundary where necessary.
 * \param[in] arr      source array (grid or data)
 * \param[in] size     size of the array
 * \param[in] start    first index of the window (may be < 0 or overrun on the right)
 * \param[in] count    number of points in the window
 * \param[in,out] dst  destination array of size `count`
 * \param[in] on_left  true if the boundary handling is on the left side
 */
template <typename T>
void fill_window(const T* arr, int size, int start, int count, T* dst, bool on_left){
   const int last = size - 1;
   const int n_extra = on_left ? std::max(0, -start)
                               : std::max(0, (start + count - 1) - last);
   const int n_valid = count - n_extra;
   const int valid_start   = on_left ? 0 : start;
   const int dst_valid_off = on_left ? n_extra : 0;
   // copy the "interior" values, return if no extrapolation is required
   copy_entries(arr, valid_start, n_valid, dst, dst_valid_off);
   if (n_extra == 0) return;
   // do extrapolation if necessary
   const T diff   = on_left ? -(arr[1] - arr[0]) : (arr[last] - arr[last - 1]);
   const T anchor = on_left ? arr[0] : arr[last];
   for (int i = 0; i < n_extra; ++i) {
      const int step = i + 1;
      const int pos  = on_left ? (n_extra - 1 - i) : (count - n_extra + i);
      dst[pos] = anchor + diff * step; // linear extrapolation
   }
}




/**
 * \brief Extracts grid and data around evaluation point, possibly extrapolated at the boundaries.
 * \param[in] derivatives select nonzero derivatives
 * \param[in] grid        array containing the grid points
 * \param[in] data        array containing the data points
 * \param[in] gridsize    size of grid and data arrays
 * \param[in] idx         index/interval of current evaluation point (preparation site)
 * \param[in,out] g       grid points around evaluation point:  g = grid(idx-3:idx+3) <-- 7 points
 * \param[in,out] z       data points around evaluation point:  z = data(idx-2:idx+2) <-- 5 points
 * \see setup_interp_param() for layout of g and z
 */
template <typename T>
inline
void setup_g_z(int derivatives, const T* grid, const T* data, int gridsize, int idx, T* g, T* z){
   // Special case for zero derivatives: only required data from local interval
   if (derivatives == 0){
      memset(g, 0, 7 * sizeof(T));   g[3] = grid[idx];   g[4] = grid[idx+1];
      memset(z, 0, 5 * sizeof(T));   z[2] = data[idx];   z[3] = data[idx+1];
      return;
   }
   // Enough data on both sides
   if (idx >= 3 && idx <= gridsize - 4) {
      copy_entries(grid, idx - 3, 7, g, 0);
      copy_entries(data, idx - 2, 5, z, 0);
      return;
   }
   // Possibly fill up by linear extrapolation on both sides
   const bool on_left = idx < 3;
   LSPI_VERBOSE_1D( displayText("Boundary: %s side, idx = %d\n", (on_left ? "left" : "right"), idx); )
   fill_window(grid, gridsize, idx - 3, 7, g, on_left);
   fill_window(data, gridsize, idx - 2, 5, z, on_left);
}





/**
 * \brief Compute and store interpolation coefficients
 * \param[in]  derivatives  if set to 1, compute "p"; if set to 2, compute "p" and "P"
 * \param[in]  g            selection of grid points in current interpolation interval
 * \param[in]  z            selection of data points in current interpolation interval
 * \param[in]  gamma        smoothing value in interval [0,1]
 * \param[out] interp_param target array to write interpolation parameters
 * \param      dummy        added for memory alignment -- UNUSED -- approx. 50% speed up when compiled with -O3
 */
template <typename T>
void setup_interp_param(int derivatives, const T* g, const T* z, const T& gamma, T* interp_param, int dummy = 0) {
   LSPI_VERBOSE_1D( displayText("Seting up interpolation parameters for derivatives %d\n", derivatives); )
   // Data layout:
   // idx h           0                1                2              3              4               5           
   //       xi-3 <==h_i-3==> xi-2 <==h_i-2==> xi_1 <==h_i-1==> xi <==h_i==> xi+1 <==h_i+1==> xi+2 <==h_i+2==> xi+3
   // idx g  0                1                2                3            4                5                6  
   //                         |                |                |            |                |                   
   // idx z                   0                1                2            3                4                   
   //                        zi-2 <==y_i-2==> zi-1 <==y_i-2==> zi <==y_i==> zi+1 <==y_i+1==> zi+2                 
   // idx y,delta,hdm,hdp              0       |        1       |      2     |        3        
   //                              Delta_i-2   |   Delta_i-1    |   Delta_i  |   Delta_i+1     
   // hdm                         h_i-2/h_i-3  |   h_i-1/h_i-2  |  h_i/h_i-1 |   h_i+1/h_i     
   // hdp                         h_i-2/h_i-1  |   h_i-1/h_i    |  h_i/h_i+1 |  h_i+1/h_i+2    
   //                                        p_im             p_i          p_ip     <-- for derivative >= 1
   //                                        P_im             P_i          P_ip     <-- for derivative >= 2
   // INPUTS: 
   //  - g: vector of grid points (size is 7)
   //  - z: vector of data points (size is 5)
   // OUTPUTS: interp_param
   //  - h:        vector of x_i+1 - x_i
   //  - sign:     sign of z_i+1 - z_i
   //  - delta:    (z_i+1 - z_i) / h_i
   //  - hdm:      h_i / h_i-1
   //  - hdp:      h_i / h_i+1
   (void) dummy;  // mark as unused: Speed up when compiled with -O3 (compiler magic, should be analyzed somewhen)
   //===== QUICK SETUP FOR ZERO DERIVATIVES
   if (derivatives == 0) {
      interp_param[0] = g[3];      // grid[idx];
      interp_param[1] = g[4]-g[3]; // h[idx] = grid[idx+1] - grid[idx];
      interp_param[2] = z[2];      // data[idx];
      interp_param[3] = z[3];      // data[idx+1];
      interp_param[4] = (T) 0.0;
      interp_param[5] = (T) 0.0;
      interp_param[6] = (T) 0.0;
      interp_param[7] = (T) 0.0;
      return;
   }

   //===== COMPUTE DIFFERENCES 
   T h[6];  diff(g,6,h);   // grid diff
   T y[4];  diff(z,4,y);   // data diff
   LSPI_VERBOSE_1D( displayVec(h, 6, "h"); displayVec(y, 4, "y") );
   T* h_m = h;             // h_i-1
   T* h_i = h + 1;         // h_i  // notation shift
   T* h_p = h + 2;         // h_i+1
   
   //===== COMPUTE POINTWISE QUOTIENTS
   T hdm[4], hdp[4], deltas[4];
   cpq_array(h_i, h_m, 4, hdm);    // pointwise quotients
   cpq_array(h_i, h_p, 4, hdp);    //
   cpq_array(y  , h_i, 4, deltas); // delta = y ./ h
   
   //===== COMPUTE WEIGHTED HARMONIC MEANS 
   T p_im,p_i,p_ip;
   T left = 5.0/gamma - 2.0;
   //                                    delta[i-1]  delta[i]  h[i] h[i+1] hdm[i] hdm[i+1] hdp[i] hdp[i+1],gamma, left
   p_im = compute_weighted_harmonic_mean(deltas[0], deltas[1], h[1], h[2], hdm[0], hdm[1], hdp[0], hdp[1], gamma, left);
   p_i  = compute_weighted_harmonic_mean(deltas[1], deltas[2], h[2], h[3], hdm[1], hdm[2], hdp[1], hdp[2], gamma, left);
   p_ip = compute_weighted_harmonic_mean(deltas[2], deltas[3], h[3], h[4], hdm[2], hdm[3], hdp[2], hdp[3], gamma, left);
   LSPI_VERBOSE_1D( displayText( "Values:  p_im = %g   p_i = %g   p_ip = %g\n", p_im, p_i, p_ip); )
   
   //===== GET SIGNUM OF DATA DIFFERENCES (ASCENDING OR DESCENDING?)
   T signs[4];
   signs[0] = get_signum(y[0]);
   signs[1] = get_signum(y[1]);
   signs[2] = get_signum(y[2]);
   signs[3] = get_signum(y[3]);
   
   //===== CALCULATE P
   T P_im = 0, P_i = 0, P_ip = 0; // initialization to 0 needed (see below)
   if (derivatives == 2) {
      const double lintol = 1e-4; // to incorporate linear data, second derivative is set to zero at grid points.
      P_im = calc_P(p_im, h[1], h[2], signs[0], signs[1], gamma);
      P_i  = (fabs(p_i-p_im) < lintol)  ?  0.0  :  calc_P(p_i , h[2], h[3], signs[1], signs[2], gamma);
      P_ip = (fabs(p_ip-p_i) < lintol)  ?  0.0  :  calc_P(p_ip, h[3], h[4], signs[2], signs[3], gamma);
      LSPI_VERBOSE_1D( displayText("Values:  P_im = %g  P_i = %g  P_ip 0 %g\n", P_im, P_i, P_ip); )
   }

   //===== CHECK SUFFICIENT CONDITIONS FOR NONZERO DERIVATIVES
   bool invalid_im, invalid_i;
   if (derivatives == 1){
      //== FOR 1st DERIVATIVE
      invalid_im = ( signs[1]*(60.0*y[1]-24.0*h[2]*(p_im + p_i )) < 0 );
      invalid_i  = ( signs[2]*(60.0*y[2]-24.0*h[3]*(p_i  + p_ip)) < 0 );
      LSPI_VERBOSE_1D( displayText("Sufficient condition violation? -- im: %d  i: %d\n ", invalid_im, invalid_i); )
      if (invalid_im) {p_i  = 0.0;} //rule in interval i-1
      if (invalid_i ) {p_ip = 0.0;} //rule in interval i
   }else if (derivatives == 2) {
      invalid_im = ( signs[1]*(60.0*y[1]-24.0*h[2]*(p_im + p_i)+3*h[2]*h[2]*(P_im - P_i)) < 0 );
      invalid_i =  ( signs[2]*(60.0*y[2]-24.0*h[3]*(p_i  +p_ip)+3*h[3]*h[3]*(P_i  -P_ip)) < 0 );
      LSPI_VERBOSE_1D( displayText("Sufficient condition violation? -- im: %d  i: %d\n ", invalid_im, invalid_i); )
      if (invalid_im){ p_i  = 0.0; P_i  = 0.0;} //rule in interval i-1
      if (invalid_i ){ p_ip = 0.0; P_ip = 0.0;} //rule in interval i
   }

   //===== SET THE INTERPOLATION PARAMETERS
   interp_param[0] = g[3];
   interp_param[1] = h[3];   // x_i+1 - x_i
   interp_param[2] = z[2];   // z_i
   interp_param[3] = z[3];   // z_i+1
   interp_param[4] = p_i;
   interp_param[5] = p_ip;
   interp_param[6] = P_i;
   interp_param[7] = P_ip;

}



/**
 * \brief Interface to evaluate interpolation at specified point.
 * \param p Query point for interpolation
 * \param interp_param vector of interpolation coefficients [x_i,h_i,z_i,z_p1,p_i,p_ip,P_i,P_ip]
 * \param derivative Derivative to be evaluated [0,1,2]
 * \return Value of evaluated interpolation at point p
 */
template <typename T>
inline
T eval_interpolation_1D(const T& p, const T* interp_param, int derivatives){
   LSPI_VERBOSE_1D( displayText("Starting evaluation of 1D interpolation at point %g\n", p) );
   T result;
   T u = (p-interp_param[0]) / interp_param[1];  //  u = (p-g)/h  with  h = (x_i+1 - x_i)
   switch (derivatives) {
      case 1:  result = eval_weighted_hermite_basis_functions_DERIV1(u,interp_param); break;
      case 2:  result = eval_weighted_hermite_basis_functions_DERIV2(u,interp_param); break;
      case 0:  result = eval_weighted_hermite_basis_functions_DERIV0(u,interp_param); break;
      default: result = eval_weighted_hermite_basis_functions_DERIV1(u,interp_param); break; // default 1
   }
   LSPI_VERBOSE_1D( displayText("EVAL_INTERPOLATION:  p = %g  u = %g  result = %g\n", p, u, result); )
   return result;
}


/**
 * \brief Performs interpolation at requested point using quintic hermite polynomials.
 * \param[in] p           point of evaluation
 * \param[in] grid        grid of data
 * \param[in] data        data to be interpolated
 * \param[in] gridsize    size of data/grid arrays
 * \param[in] derivatives specify to use 0,1,2 derivatives
 * \param[in] gamma       smoothing parameter, within [0,1]
 * \param[in,out] interp_param interpolation coefficients - reusable storage/cache
 * \param[in] idx_old     previously queried index/grid interval - belongs to interp_param
 * \return Interpolated data at point p
 */
template <typename T>
T compute_interpolation_1D(const T& p, const T* grid, const T* data, int gridsize, int derivatives, const T& gamma, T* interp_param, int& idx_old){
   LSPI_VERBOSE_1D( displayText("Starting index search ... "); )
   int idx = find_index_from(p, grid, gridsize, idx_old); // search interval in grid that contains p
   if (idx_old != idx) {  // index/interval changed: interpolation coefficients must be updated
      idx_old = idx;      // store idx
      LSPI_VERBOSE_1D( displayText("Index has changed. Updating interpolation parameters ... "); )
      T g[7]{}, z[5]{};
      setup_g_z(derivatives, grid, data, gridsize, idx, g, z);
      setup_interp_param(derivatives, g, z, gamma, interp_param);
   }
   T result = eval_interpolation_1D(p, interp_param, derivatives);
   LSPI_VERBOSE_1D( displayText("Interpolation result: %g\n", result); )
   return result;
}


} // end of namespace


#endif // include guard
