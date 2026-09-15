// THIS FILE CONTAINS OLD LEFTOVERS

template <typename T>
inline
T eval_interpolation_DERIV1(const T& p, const T*interp_param){
   // function for evalution of interpolation
   // INPUT:
   //  -p: point of evaluation
   //  -interp_param = [x,h_i,z[2],z[3],p_i,p_ip,P_i,P_ip] 
   // OUTPUT:
   //  -value of evaluted interpolation at point p
   T u = (p-interp_param[0]) / interp_param[1];
   T result = eval_weighted_hermite_basis_functions_DERIV1(u,interp_param);
   LSPI_VERBOSE_1D(std::cout<<"INTERP_DERIV1: p = "<<p<<" : value = "<<result<<std::endl;)
   return result;
} 



template <typename T>
inline
T eval_interpolation_DERIV2(const T& p,const T*  interp_param){
   // function for evalution of interpolation
   // INPUT:
   //  -p: point of evaluation
   //  -interp_param = [x,h_i,z[2],z[3],p_i,p_ip,P_i,P_ip] 
   // OUTPUT:
   //  -value of evaluted interpolation at point p
   T u = (p-interp_param[0]) / interp_param[1];
   T result = eval_weighted_hermite_basis_functions_DERIV2(u, interp_param);
   LSPI_VERBOSE_1D(std::cout<<"INTERP_DERIV2: p = "<<p<<" : value = "<<result<<std::endl;)
   return result;
} 



/**
 * \brief Compute weighted harmonic mean (whm)
 */
template <typename T>
T compute_weighted_harmonic_mean(const T& delta_im, const T& delta_i,
                                 const T& h_im    , const T& h_i,
                                 const T& hdm_im  , const T& hdm_i,
                                 const T& hdp_im  , const T& hdp_i,
                                 const T& alpha   , const T& left)
{
   // function to decide, if at position i we have to take the whm or 0.
   // conditions are as follows:
   // ensure not to divide by zero
   LSPI_VERBOSE_1D(std::cout<< "info for whm (hdm_i,hdp_im): "<< hdm_i<< " "<<hdp_im;)
   bool c1 = (delta_im*delta_i>0);
   // checks on the grid to get a possible solution
   bool c2 = (left>=hdm_im+hdp_im);
   bool c3 = (left>=hdm_i+hdp_i);
   // checks, that we satisfy the cufficient condition for monoticity preserving,
   // when the rule in STEP 3 in function eval_interpolation_DERIV1 is applied.
   bool c4 = (left>=2*hdm_i);  // 3 - 2 * h_i/h_i-1 >=0
   bool c5 = (left>=2*hdp_im); // 3 - 2 * h_i-1/h_i >=0
   LSPI_VERBOSE_1D(std::cout << "conditions for whm: "<< c1 << " "<< c2 << " "<< c3 << " "<< c4 << " "<< c5 << std::endl;)
   if(c1 && c2 && c3 && c4 && c5){
      //Formula:   whm = (h_i-1 + h_i) / (h_i/delta_i-1 + h_i-1/delta_i)
      return (T) alpha * (h_im  + h_i) / (h_i/delta_im  + h_im /delta_i);
   }else{
      return (T) 0.0;
   }
}




/**
 * \brief Compute weighted harmonic means (whm) or decide to use 0.
 * \input[in]
 */
template <typename T>
T compute_weighted_harmonic_means(const T& deltas, const T& h, const T& hdm, const T& hdp, const T& alpha, const T& left,
                                  T& p_im, T& p_i, T& p_ip)
{
    // I don't think this can be optimized better than the compiler does...
}




/**
 * \brief Compute and store interpolation coefficients
 * \param[in]  p point of evaluation/interpolation
 * \param[in]  x array of grid points
 * \param[in]  z array of data points
 * \param[in]  alpha value in interval [0,1]
 * \param[out] interp_param target array to write interpolation coefficients
 * \param[in]  derivative if set to 1, compute "p"; if set to 2, compute "p" and "P"
 */
template <typename T>
void setup_interpolation_coeffs(const T* x, const T* z, const T& alpha, T* interp_param, int derivatives) {
   LSPI_VERBOSE_1D(std::cout<<"Evaluation of interpolation at point "<<p<<" has started"<<std::endl;)
   // Data layout:
   // idx h           0                1                2              3              4               5           
   //       xi-3 <==h_i-3==> xi-2 <==h_i-2==> xi_1 <==h_i-1==> xi <==h_i==> xi+1 <==h_i+1==> xi+2 <==h_i+2==> xi+3
   // idx x  0                1                2                3            4                5                6  
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
   //  - x: vector of grid points (size is 7)
   //  - z: vector of data points (size is 5)
   // OUTPUTS: interp_param
   //  - h:        vector of x_i+1 - x_i
   //  - sign:     sign of z_i+1 - z_i
   //  - delta:    (z_i+1 - z_i) / h_i
   //  - hdm:      h_i / h_i-1
   //  - hdp:      h_i / h_i+1
   
   //===== COMPUTE DIFFERENCES 
   T h[6];  diff(x,6,h);   // grid diff
   T y[4];  diff(z,4,y);   // data diff
   LSPI_VERBOSE_1D(std::cout<< "h vector with "; print_vector(h,6); std::cout<< "y vector with "; print_vector(y,4);)
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
   T left = 5.0/alpha -2.0;
   //                                    delta[i-1]  delta[i]  h[i] h[i+1] hdm[i] hdm[i+1] hdp[i] hdp[i+1],alpha, left
   p_im = compute_weighted_harmonic_mean(deltas[0], deltas[1], h[1], h[2], hdm[0], hdm[1], hdp[0], hdp[1], alpha, left);
   p_i  = compute_weighted_harmonic_mean(deltas[1], deltas[2], h[2], h[3], hdm[1], hdm[2], hdp[1], hdp[2], alpha, left);
   p_ip = compute_weighted_harmonic_mean(deltas[2], deltas[3], h[3], h[4], hdm[2], hdm[3], hdp[2], hdp[3], alpha, left);
   LSPI_VERBOSE_1D(std::cout<<"Evaluation of p_im, p_i, p_ip is "<<p_im<< " "<< p_i << " "<<p_ip<< std::endl;)
   
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
      P_im = calc_P(p_im, h[1], h[2], signs[0], signs[1], alpha);
      P_i  = (fabs(p_i-p_im) < lintol)  ?  0.0  :  calc_P(p_i , h[2], h[3], signs[1], signs[2], alpha);
      P_ip = (fabs(p_ip-p_i) < lintol)  ?  0.0  :  calc_P(p_ip, h[3], h[4], signs[2], signs[3], alpha);
      LSPI_VERBOSE_1D(std::cout<<"P_im = "<<P_im<<" P_i = "<<P_i<< " P_ip = "<<P_ip<< std::endl;)
   }

   //===== CHECK SUFFICIENT CONDITIONS FOR NONZERO DERIVATIVES
   bool invalid_im, invalid_i;
   if (derivatives == 1){
      //== FOR 1st DERIVATIVE
      invalid_im = ( signs[1]*(60.0*y[1]-24.0*h[2]*(p_im + p_i )) < 0 );
      invalid_i  = ( signs[2]*(60.0*y[2]-24.0*h[3]*(p_i  + p_ip)) < 0 );
      LSPI_VERBOSE_1D(std::cout<<"Sufficient condition violation? -- im: "<<invalid_im<<" i:"<<invalid_i<<std::endl;)
      if (invalid_im) {p_i  = 0.0;} //rule in interval i-1
      if (invalid_i ) {p_ip = 0.0;} //rule in interval i
   }else if (derivatives == 2) {
      invalid_im = ( signs[1]*(60.0*y[1]-24.0*h[2]*(p_im + p_i)+3*h[2]*h[2]*(P_im - P_i)) < 0 );
      invalid_i =  ( signs[2]*(60.0*y[2]-24.0*h[3]*(p_i  +p_ip)+3*h[3]*h[3]*(P_i  -P_ip)) < 0 );
      LSPI_VERBOSE_1D(std::cout<<"Sufficient condition violation  -- im: "<<invalid_im<<" i:"<<invalid_i<<std::endl;)
      if (invalid_im){ p_i  = 0.0; P_i  = 0.0;} //rule in interval i-1
      if (invalid_i ){ p_ip = 0.0; P_ip = 0.0;} //rule in interval i
   }

   //===== SET THE INTERPOLATION COEFFICIENTS
   interp_param[0] = x[3];
   interp_param[1] = h[3];   //
   interp_param[2] = z[2];   // z_i
   interp_param[3] = z[3];   // z_i+1
   interp_param[4] = p_i;
   interp_param[5] = p_ip;
   interp_param[6] = P_i;
   interp_param[7] = P_ip;
}



template <typename T>
inline
T eval_interpolation_DERIV0(const T& p, const T& x0, const T& x1, const T& z0, const T& z1){
   // function for simplified evalution of interpolation
   // INPUT:
   //  -p: point of evaluation
   //  -x: vector of grid points (x_i,x_i+1)
   //  -z: vector of data points (z_i,z_i+1)
   // OUTPUT:
   //  -value of evaluted interpolation at point p
   LSPI_VERBOSE_1D(std::cout<<"Evaluation of interpolation at point "<<p<<" has started"<<std::endl;)
   T h = x1-x0;
   T u = (p-x0)/h;
   T result = eval_weighted_hermite_basis_functions_DERIV0(u,z0,z1);
   LSPI_VERBOSE_1D( std::cout<<"INTERP_DERIV0:  h = "<<h<<" -- u = "<<u<<" -- p = "<<p<<" : value = "<<result<<std::endl; )
   return result;
}





/**
 * \brief Evaluate weighted hermite basis functions -- with NO DERIVATIVE weights
 * \param[in] u   point of evaluation in the respective interval (p-x_i)/h_i
 * \param[in] z0  parameter z_i
 * \param[in] z1  parameter z_i+1
 */
template <typename T>
inline
T eval_weighted_hermite_basis_functions_DERIV0(const T& u, const T& z0, const T& z1){
   T mu  = 1-u;
   T u2  = u*u;
   T umu = u*mu;
   T row = umu *(1+u-9*u2+6*u2*u);
   return (T) (z0*(mu + row) + z1*(u - row));
}






/**
 * \brief Extracts grid and data around evaluation point, possibly extrapolated at the boundaries.
 * \param[in] grid      array containing the grid points
 * \param[in] data      array containing the data points
 * \param[in] gridsize  size of grid and data arrays
 * \param[in] idx       index/interval of current evaluation point (preparation site)
 * \param[in,out] x     grid points around evaluation point:  x = grid(idx-3:idx+3) <-- 7 points
 * \param[in,out] z     data points around evaluation point:  z = data(idx-2:idx+2) <-- 5 points
 */
template <typename T>
void setup_g_z(const T* grid, const T* data, int gridsize, int idx, T* x, T* z){
   if (idx>=3 && idx<=gridsize-4) {  // we have enough data, i.e. in the middle of grid and data
      copy_entries(grid, idx-3 , 7, x, 0);
      copy_entries(data, idx-2 , 5, z, 0);
   }else if(idx < 3){ // we have to add data points on the left side
      int m_idx = 3 - idx;
      LSPI_VERBOSE_1D(std::cout<<"Boundary on the left side with m_idx="<<m_idx<<std::endl;)
      copy_entries(grid, 0, 7-m_idx, x, m_idx);
      T x_diff = grid[1] - grid[0];
      for(int idx_int=1; idx_int<=m_idx; ++idx_int){
         x[m_idx - idx_int] = grid[0] - x_diff*idx_int; // linear extrapolation on the left side of data
      }
      if(m_idx == 0){
         copy_entries(data, 0, 5, z, 0);
      }else{
         copy_entries(data, 0, 5-m_idx+1, z, m_idx-1);  // same for data if needed
         T z_diff = data[1] - data[0]; 
         for(int idx_int=0; idx_int<m_idx-1; ++idx_int){
            z[m_idx-2 - idx_int] = data[0] - z_diff*(idx_int+1);
         }   
      }
   }else{  // we have to add data points on the right side
      int p_idx = idx - (gridsize-4);  
      LSPI_VERBOSE_1D(std::cout<<"Boundary on the right side with p_idx="<<p_idx<<std::endl;)
      copy_entries(grid, gridsize-7+p_idx, 7-p_idx, x, 0);
      int last = gridsize-1; 
      T x_diff = grid[last]-grid[last-1];
      for(int idx_int=0; idx_int<p_idx; ++idx_int){  
            x[7-p_idx+idx_int] = grid[last] + x_diff*(idx_int+1); // linear extrapolation on the right side
      }
      // special cases
      if(p_idx == 1){
            copy_entries(data, gridsize-5, 5, z, 0);
      }else if (p_idx == 2){
            copy_entries(data, gridsize-4, 4, z, 0);
            z[4] = data[last] + (data[last]-data[last-1]); // last value in z is missing: need to extrapolate
      }else{ 
            LSPI_VERBOSE_1D(std::cout<<"WARNING: p_idx>2. This case should not exist. I'll take last 5 values from data";)
            copy_entries(data, gridsize-6, 5, z, 0);
      }
   }
   // TODO:: I bet the extrapolation cases can be easily unified 
}




/**
 * \brief Extracts grid and data around evaluation point, possibly extrapolated at the boundaries.
 * \param[in] grid      array containing the grid points
 * \param[in] data      array containing the data points
 * \param[in] gridsize  size of grid and data arrays
 * \param[in] idx       index/interval of current evaluation point (preparation site)
 * \param[in,out] x     grid points around evaluation point:  x = grid(idx-3:idx+3) <-- 7 points
 * \param[in,out] z     data points around evaluation point:  z = data(idx-2:idx+2) <-- 5 points
 */
template <typename T>
void setup_g_z(const T* grid, const T* data, int gridsize, int idx, T* x, T* z){
   if (idx >= 3 && idx <= gridsize - 4) {  // enough data on both sides
      copy_entries(grid, idx - 3, 7, x, 0);
      copy_entries(data, idx - 2, 5, z, 0);
      return;
   }
   const bool on_left = idx < 3;
   LSPI_VERBOSE_1D(std::cout << "Boundary " << (on_left ? "left" : "right") << " side, idx=" << idx << std::endl;)
   fill_window(grid, gridsize, idx - 3, 7, x, on_left);
   fill_window(data, gridsize, idx - 2, 5, z, on_left);
}



enum Direction {
    X, Y
};


/**
 * \brief 
 */
template <typename T>
inline
void copy_entries_from_two_dim(Direction dir, const T* data, int idxRunning, int count, int idxFixed, int size_x, T* ar2, T* ar3, int idx_begin_out, bool x_direction ){
    //function to copy from data into tow arrays.
    // IN:
    // - data with data(i,j) = data[i*size_x+j], i.e. with fixed i we get in x direction data and with fixed j we get the data along the y direction
    // - idx_begin_ar1
    if(dir == Direction::X){
        for(int j = 0; j<count; ++j){
            ar2[idx_begin_out+j] = LSPI_MATRIX_ELEMENT(data, idxFixed  , idxRunning+j, size_x);
            ar3[idx_begin_out+j] = LSPI_MATRIX_ELEMENT(data, idxFixed+1, idxRunning+j, size_x);
        }
    }else{
        for(int i = 0; i<count; ++i){
            ar2[idx_begin_out+i] = LSPI_MATRIX_ELEMENT(data, idxRunning+i, idxFixed  , size_x);
            ar3[idx_begin_out+i] = LSPI_MATRIX_ELEMENT(data, idxRunning+i, idxFixed+1, size_x);
        }    
    }
}






template <typename T>
void setup_g_zg_zgh(const T* grid, const T* data, int gridsize, int leadim, int idx, int idxFixed, 
                           T* g, T * z_g, T* z_gh, bool x_direction){
 

    if (idx>=3 && idx<=gridsize-4){ 
        // we have enough data to the left and to the right
        copy_entries(grid, idx-3, 7, g);
        copy_entries_from_two_dim(data, idx-2, 5, idxFixed, leadim, z_g, z_gh, 0, x_direction);
    }else if(idx <3){ 
        // we have to add data points on the left side
        int m_idx = 3-idx;
        LSPI_VERBOSE_2D(std::cout<<"Boundary on the left side with m_idx="<<m_idx<<std::endl;)
        copy_entries(grid, 0, 7-m_idx, g, m_idx);
        T x_diff = grid[1]-grid[0];
        for(int idx_int = 1; idx_int<=m_idx;idx_int++){
            g[m_idx-idx_int] = grid[0]-x_diff*idx_int; //linear extrapolation
        }
        if(m_idx == 0){
            copy_entries_from_two_dim(data, 0, 5, idxFixed, leadim, z_g, z_gh, 0      , x_direction);
        }else{
            copy_entries_from_two_dim(data, 0, 5, idxFixed, leadim, z_g, z_gh, m_idx-1, x_direction);
            if(x_direction){
                T z_diff   = LSPI_MATRIX_ELEMENT(data, idxFixed  , 1, leadim) - LSPI_MATRIX_ELEMENT(data, idxFixed  , 0, leadim);  //for x_direction for fixed row i
                T z_p_diff = LSPI_MATRIX_ELEMENT(data, idxFixed+1, 1, leadim) - LSPI_MATRIX_ELEMENT(data, idxFixed+1, 0, leadim);
                for(int idx_int = 0; idx_int<m_idx-1;idx_int++){
                    z_g [m_idx-2 - idx_int] = LSPI_MATRIX_ELEMENT(data, idxFixed  , 0, leadim) -   z_diff*(idx_int+1); //write missing values
                    z_gh[m_idx-2 - idx_int] = LSPI_MATRIX_ELEMENT(data, idxFixed+1, 0, leadim) - z_p_diff*(idx_int+1);
                }
            }else{
                T z_diff   = LSPI_MATRIX_ELEMENT(data, 1, idxFixed  , leadim) - LSPI_MATRIX_ELEMENT(data, 0, idxFixed  ,leadim); //for y_direction for fixed column j
                T z_p_diff = LSPI_MATRIX_ELEMENT(data, 1, idxFixed+1, leadim) - LSPI_MATRIX_ELEMENT(data, 0, idxFixed+1,leadim);
                for(int idx_int = 0; idx_int<m_idx-1;idx_int++){  
                    z_g [m_idx-2 - idx_int] = LSPI_MATRIX_ELEMENT(data, 0, idxFixed  , leadim)   -z_diff*(idx_int+1); //write missing values
                    z_gh[m_idx-2 - idx_int] = LSPI_MATRIX_ELEMENT(data, 0, idxFixed+1, leadim) -z_p_diff*(idx_int+1);
                }
            }

        }
    }else{
        // we have to add data points on the right side
        int p_idx = idx-(gridsize-4);
        LSPI_VERBOSE_2D(std::cout<<"Boundary on the right side with p_idx="<<p_idx<<std::endl;)
        copy_entries(grid,gridsize-7+p_idx,7-p_idx,g);
        int last = gridsize-1;
        T x_diff = grid[last]-grid[last-1];
        for(int idx_int = 0; idx_int<p_idx; ++idx_int){  
            g[7-p_idx+idx_int] = grid[last] + x_diff*(idx_int+1); //linear extrapolation
        }
        if(p_idx == 1){
            copy_entries_from_two_dim(data,gridsize-5,5,idxFixed,leadim,z_g,z_gh,0,x_direction);
        }else if (p_idx == 2){    //if needed also with z and zh
            copy_entries_from_two_dim(data,gridsize-4,4,idxFixed,leadim,z_g,z_gh,0,x_direction); // last value in z is missing
            if(x_direction){
                z_g[4]   = LSPI_MATRIX_ELEMENT(data, idxFixed  , last, leadim) + (LSPI_MATRIX_ELEMENT(data, idxFixed  , last, leadim) - LSPI_MATRIX_ELEMENT(data, idxFixed  , (last-1), leadim)); //for fixed  row i
                z_gh[4]  = LSPI_MATRIX_ELEMENT(data, idxFixed+1, last, leadim) + (LSPI_MATRIX_ELEMENT(data, idxFixed+1, last, leadim) - LSPI_MATRIX_ELEMENT(data, idxFixed+1, (last-1), leadim));
            }else{
                z_g[4]   = LSPI_MATRIX_ELEMENT(data, last, idxFixed  , leadim) + (LSPI_MATRIX_ELEMENT(data, last, idxFixed  , leadim) - LSPI_MATRIX_ELEMENT(data, last-1, idxFixed  ,leadim)); //for fixed column j
                z_gh[4]  = LSPI_MATRIX_ELEMENT(data, last, idxFixed+1, leadim) + (LSPI_MATRIX_ELEMENT(data, last, idxFixed+1, leadim) - LSPI_MATRIX_ELEMENT(data, last-1, idxFixed+1,leadim));
            }
        }else{ 
            LSPI_VERBOSE_2D( std::cout<<"WARNING: p_idx>2. This case should not exist, take last 5 values from data"; );
            copy_entries_from_two_dim(data,gridsize-6,5,idxFixed,leadim,z_g,z_gh,0,x_direction);
        }
    }
}






#define LSPI_COPY_X(source, start, count, target )

/**
 * \brief Extracts several matrix elements from specified row and next row into target arrays.
 * \param[in] data Data matrix (data to be interpolated)
 * \param[in] idx_begin
 */
template <typename T>
inline
void copy_entries_2D_x( const T* data, int idx_begin, int count, int idxFixed, int leadim, T* ar2, T* ar3, int idx_begin_out)
{
    for(int i = 0; i<count; ++i){
        ar2[idx_begin_out+i] = LSPI_MATRIX_ELEMENT(data, idxFixed  , idx_begin+i, leadim); // for x_dir for fixed row i
        ar3[idx_begin_out+i] = LSPI_MATRIX_ELEMENT(data, idxFixed+1, idx_begin+i, leadim); // for x_dir for fixed row i+1
    }
}



/**
* \brief Extracts several matrix elements from specified column and next column into target arrays.
 */
template <typename T>
inline
void copy_entries_2D_y( const T* data, int idx_begin, int count, int idxFixed, int leadim, T* ar2, T* ar3, int idx_begin_out, bool x_direction ){
    //function to copy from data into tow arrays.
    // IN:
    // - data with data(i,j) = data[i*leadim+j], i.e. with fixed i we get in g direction data and with fixed j we get the data along the y direction
    // - idx_begin_ar1
    for(int i = 0; i<count; ++i){
        ar2[idx_begin_out+i] = LSPI_MATRIX_ELEMENT(data, idx_begin+i, idxFixed  , leadim); //for y_dir for fixed column j
        ar3[idx_begin_out+i] = LSPI_MATRIX_ELEMENT(data, idx_begin+i, idxFixed+1, leadim); //for y_dir for fixed column j+1
    }    
}