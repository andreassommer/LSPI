#ifndef LSPI_2D_HEADER
#define LSPI_2D_HEADER

/**
 *  \file Header of 2D local shape preserving interpolation 
 *  \author Julian Niederer <julian.niederer@iwr.uni-heidelberg.de>
 *  \author Andreas Sommer <code@andreas-sommer.eu>
 */



// 2D version relies on 1D local interpolation
#include "LSPI_1D.hpp"



namespace LSPI {



/**
 * \brief Evaluate 2D interpolation formula at specified point.
 * \param[in] patchdata data for interpolation: [g, y, g+h_x, y+h_y, z_ll, z_ul, z_lr, z_ur]
 * \param[in] current   vector: [x_cur, y_cur, lower, upper, left, right] (built by 1D interpolation)
 * \return Value of 2D interpolation at point (x_cur,y_cur)
 */
template <typename T>
T eval_interpolation_2D(const T* patchdata, const T* current)
{
    // - patchdata:  [g, y, g+h_x, y+h_y, z_ll, z_ul, z_lr, z_ur]
    //           0  1  2      3      4     5     6     7
    //   --- changes only if query points leaves the previous patch
    //
    // - current:  [x_cur, y_cur, lower, upper, left, right]
    //              0      1      2      3      4     5
    //   --- changes for every query point
    //
    //  current: 
    //  0: x_cur:  query point (g direction)
    //  1: y_cur:  query point (y direction)
    //  2: lower:  data value at (x_cur, y    ), 1d-interpolated
    //  3: upper:  data value at (x_cur, y+h_y), 1d-interpolated
    //  4: left :  data value at (g    , y_cur), 1d-interpolated
    //  5: right:  data value at (g+h_x, y_cur), 1d-interpolated
    //
    //                               (x_cur,y+h_y)
    //         (g,y+h_y) + ------------*--------------------------+ (g+h_x,y+h_y) 
    //                   | z_ul      upper                   z_ur |
    //                   |             :                          |
    //         (g,y_cur) * left .....  g <-(x_cur,y_cur)... right * (g+h_x,y_cur)
    //                   |             :                          |
    //                   | z_ll      lower                   z_lr |
    //             (g,y) +-------------*--------------------------+ (g+h_x,y)
    //                               (x_cur,y)
    //
    T u_x = (current[0]-patchdata[0])/(patchdata[2]-patchdata[0]);  // scale into
    T u_y = (current[1]-patchdata[1])/(patchdata[3]-patchdata[1]);  // interval [0,1]
    T weight_x = smoothstep5(u_x);
    T weight_y = smoothstep5(u_y);
    LSPI_VERBOSE_2D( displayText("weight_x = %g   weight_y = %g\n", weight_x, weight_y); )
    return (T) (weight_y    *   weight_x   *(current[2] + current[4] - patchdata[4])    // (lower + left  + z_ll)
            + (1.0-weight_y)*   weight_x   *(current[3] + current[4] - patchdata[5])    // (upper + left  + z_ul)
            +   weight_y    *(1.0-weight_x)*(current[2] + current[5] - patchdata[6])    // (lower + right + z_lr)
            + (1.0-weight_y)*(1.0-weight_x)*(current[3] + current[5] - patchdata[7]) ); // (upper + right + z_ur)
}




/**
 * \brief Extracts several matrix elements from specified row and next row into target arrays.
 * \param[in] dir        Direction X or Y
 * \param[in] data       Data matrix (data to be interpolated)
 * \param[in] idxRunning The "running" index start index (for travering rows or cols), 
 * \param[in] count      Number of elements to be copied
 * \param[in] idxFixed   The "fixed" index start index (to copy row and row+1  or  col and col+1)
 * \param[in] leadim     leading dimension for flat memory (data matrix, row-major order, equals gridXsize)
 * \param[out] data0     array containing the extracted data for specified fixed row/col
 * \param[out] data1     array containing the extracted data for specified fixed row/col + 1
 * \param[out] idxOut    index in data0 and data1 to start writing outputs
 */
template <typename T>
inline
void copy_entries_2D(Direction dir, const T* data, int idxRunning, int count, int idxFixed, int leadim, T* ar2, T* ar3, int idxOut){
    // Note: data(i,j) = data[i*leadim+j] -- fix i for x direction, fix j for y direction
    if(dir == Direction::X){
        for(int j = 0; j<count; ++j){
            ar2[idxOut+j] = LSPI_MATRIX_ELEMENT(data, idxFixed  , idxRunning+j, leadim);
            ar3[idxOut+j] = LSPI_MATRIX_ELEMENT(data, idxFixed+1, idxRunning+j, leadim);
        }
    }else{
        for(int i = 0; i<count; ++i){
            ar2[idxOut+i] = LSPI_MATRIX_ELEMENT(data, idxRunning+i, idxFixed  , leadim);
            ar3[idxOut+i] = LSPI_MATRIX_ELEMENT(data, idxRunning+i, idxFixed+1, leadim);
        }    
    }
}









/**
 * \brief Retrieve data from patch, possibly extrapolate at the borders
 * \param[in] data      full data array (matrix)
 * \param[in] grid      full grid in specified direction
 * \param[in] gridsize  size of gridX
 * \param[in] leadim    leading dimension for data array
 * \param[in] idxGrid   index of patch begin in grid
 * \param[in] idxFixed  index of fixed row/column
 * \param[out] g        grid window, extracted grid points around specified index
 * \param[out] z_g      data at extracted grid window
 * \param[out] z_gh     data at extracted grid window plus one row/column
 * \param[in] dir       direction, either X or Y (enums)
 */
template <typename T>
void setup_g_zg_zgh(Direction dir, const T* grid, const T* data, int gridsize, int leadim, int idxGrid, int idxFixed, 
                    T* g, T* z_g, T* z_gh) 
{
    // Data layout
    //                               (x_cur,y+h_y)
    //         (g,y+h_y) + ------------*--------------------------+ (g+h_x,y+h_y) 
    //                   | z_ul      upper                   z_ur |
    //                   | ______        :              _________ |
    //         (g,y_cur) * |LEFT| .....  g .............| RIGHT | * (g+h_x,y_cur)
    //                   |             :                          |
    //                   | z_ll      lower                   z_lr |
    //             (g,y) +-------------*--------------------------+ (g+h_x,y)
    //                               (x_cur,y)
                    
    // elem(offset, pos) retrieves data at (idxFixed+offset, pos) or (pos, idxFixed+offset)
    // offset 0 -> z_g,   offset 1 -> z_gh
    auto elem = [&](int offset, int pos) -> T {
        return (dir==Direction::X) ? LSPI_MATRIX_ELEMENT(data, idxFixed + offset, pos, leadim)
                                   : LSPI_MATRIX_ELEMENT(data, pos, idxFixed + offset, leadim);
    };
    T* z[2] = { z_g, z_gh };

    if (idxGrid >= 3 && idxGrid <= gridsize - 4) {   // patch lies fully inside data, no extrapolation needed
        copy_entries(grid, idxGrid - 3, 7, g, 0);
        copy_entries_2D(dir, data, idxGrid - 2, 5, idxFixed, leadim, z_g, z_gh, 0);
    } else if (idxGrid < 3) {
        int m_idx = 3 - idxGrid;
        LSPI_VERBOSE_2D( displayText("Boundary on the left side with m_idx = %d\n", m_idx); )
        copy_entries(grid, 0, 7 - m_idx, g, m_idx);  // copy grid data
        T x_diff = grid[1] - grid[0];
        for (int i = 1; i <= m_idx; ++i) {
            g[m_idx - i] = grid[0] - x_diff * i;     // extrapolate missing data
        }
        int missing = (m_idx == 0) ? 0 : m_idx - 1;
        copy_entries_2D(dir, data, 0, 5-missing, idxFixed, leadim, z_g, z_gh, missing);
        for (int k = 0; k < 2 && missing > 0; ++k) {
            T diff = elem(k, 1) - elem(k, 0);
            for (int i = 0; i < missing; ++i) {
                z[k][missing - 1 - i] = elem(k, 0) - diff * (i + 1);
            }
        }
    } else {
        int p_idx = idxGrid - (gridsize - 4);
        LSPI_VERBOSE_2D( displayText("Boundary on the right side with p_idx = %d\n", p_idx); )
        copy_entries(grid, gridsize - 7 + p_idx, 7 - p_idx, g, 0);
        int last = gridsize - 1;
        T x_diff = grid[last] - grid[last - 1];
        for (int i = 0; i < p_idx; ++i) {
            g[7 - p_idx + i] = grid[last] + x_diff * (i + 1);
        }
        if (p_idx == 1) {
            copy_entries_2D(dir, data, gridsize - 5, 5, idxFixed, leadim, z_g, z_gh, 0);
        } else if (p_idx == 2) {
            copy_entries_2D(dir, data, gridsize - 4, 4, idxFixed, leadim, z_g, z_gh, 0);
            for (int k = 0; k < 2; ++k) {
                z[k][4] = elem(k, last) + (elem(k, last) - elem(k, last - 1));
            }
        } else {
            LSPI_VERBOSE_2D( displayText("WARNING: p_idx = %d > 2. This case should not exist. I'll take last 5 values from data", p_idx); )
            copy_entries_2D(dir, data, gridsize - 6, 5, idxFixed, leadim, z_g, z_gh, 0);
        }
    }
    // The both corner cases might be unified.
}






/**
 * \brief  Performs interpolation at requested point using underlying 1D interpolation.
 * \param[in] qX          query point's X value
 * \param[in] qY          query point's Y value
 * \param[in] gridX       grid in X direction
 * \param[in] gridY       grid in Y direction
 * \param[in] data        data to be interpolated
 * \param[in] gridXsize   size gridX array (and number of columns of data) = leadim for row-major order
 * \param[in] gridYsize   size gridY array (and number of rows of data)
 * \param[in] derivatives specify to use 0,1,2 derivatives
 * \param[in] gammaX      smoothing parameter for X direction, within [0,1]
 * \param[in] gammaY      smoothing parameter for Y direction, within [0,1]
 * \param[in,out] ipparam_xz   interpolation parameters for "lower"
 * \param[in,out] ipparam_xzh  interpolation parameters for "upper"
 * \param[in,out] ipparam_yz   interpolation parameters for "left"
 * \param[in,out] ipparam_yzh  interpolation parameters for "right"
 */
template <typename T>
T compute_interpolation_2D(const T& qX, const T& qY, const T* gridX, const T* gridY, const T* data, 
                           int gridXsize, int gridYsize, int derivatives, const T gammaX, const T gammaY, 
                           T* ipparam_xz, T* ipparam_xzh, int& idx_x_old,
                           T* ipparam_yz, T* ipparam_yzh, int& idx_y_old)
{
    // Interpolate data at query point q = (x_cur,y_cur).
    //                               (x_cur,y+h_y)
    //         (g,y+h_y) + ------------*--------------------------+ (g+h_x,y+h_y) 
    //                   | z_ul      upper                   z_ur |
    //                   |             :                          |
    //         (g,y_cur) * left .....  g <-p .............. right * (g+h_x,y_cur)
    //                   |             :                          |
    //                   | z_ll      lower                   z_lr |
    //             (g,y) +-------------*--------------------------+ (g+h_x,y)
    //                               (x_cur,y)

    // INDEX SEARCH
    int idx_x = find_index_from(qX, gridX, gridXsize, idx_x_old);
    int idx_y = find_index_from(qY, gridY, gridYsize, idx_y_old);

    // UPDATE PATCHDATA IF NOT IN SAME PATCH
    if(idx_x_old != idx_x || idx_y_old !=idx_y){
        LSPI_VERBOSE_2D( displayText("Entering new patch. old idx  x:%d  y:%d  ---  Updating patchdata\n", idx_x_old, idx_y_old); )
        idx_x_old = idx_x;
        idx_y_old = idx_y;
        T x[7], z_x[5], z_xh[5], y[7], z_y[5], z_yh[5];
        setup_g_zg_zgh(Direction::X, gridX, data, gridXsize, gridXsize, idx_x, idx_y, x, z_x, z_xh);
        setup_interp_param(derivatives, x, z_x , gammaX, ipparam_xz );
        setup_interp_param(derivatives, x, z_xh, gammaX, ipparam_xzh);
        setup_g_zg_zgh(Direction::Y, gridY, data, gridYsize, gridXsize, idx_y, idx_x, y, z_y, z_yh);
        setup_interp_param(derivatives, y, z_y , gammaY, ipparam_yz );
        setup_interp_param(derivatives, y, z_yh, gammaY, ipparam_yzh);
        LSPI_VERBOSE_2D( 
            displayText("Updated patchdata:\n");
            displayVec(x   , 7, "x   ");
            displayVec(z_x , 5, "z_x ");
            displayVec(z_xh, 5, "z_xh");
            displayVec(y   , 7, "y   ");
            displayVec(z_y , 5, "z_y ");
            displayVec(z_yh, 5, "z_yh");
        )
    }

    // DEFINE LSPI_MATRIX_ELEMENT AND EVAL 1D INTERPOLATION
    T lower, upper, left, right;
    LSPI_VERBOSE_2D(
        displayFmt("idx  x:%d  y:%d\n", idx_x, idx_y);
        displayVec(ipparam_xzh ,6 , "ipparam_xzh ");
        displayVec(ipparam_xz  ,6 , "ipparam_xz  ");
        displayVec(ipparam_yz  ,6 , "ipparam_yz  ");
        displayVec(ipparam_yzh ,6 , "ipparam_yzh ");
    )
    lower = eval_interpolation_1D(qX, ipparam_xz , derivatives);
    upper = eval_interpolation_1D(qX, ipparam_xzh, derivatives);
    left  = eval_interpolation_1D(qY, ipparam_yz , derivatives);
    right = eval_interpolation_1D(qY, ipparam_yzh, derivatives);
    LSPI_VERBOSE_2D(
        displayVar(lower, "lower", "  ");
        displayVar(upper, "upper", "  ");
        displayVar(left , "left ", "  ");
        displayVar(right, "right", "\n");
    )
    // WEIGHTING TO QUERY POINT q
    // data_input: [ g, y, h_x, h_y, z_ll, z_ul, z_lr, z_ur] (built from grid and data)
    //        idx:   0  1  2    3    4     5     6     7 
    T data_input[8] = {gridX[idx_x], gridY[idx_y], gridX[idx_x+1], gridY[idx_y+1],
                       LSPI_MATRIX_ELEMENT(data, idx_y  , idx_x  , gridXsize),
                       LSPI_MATRIX_ELEMENT(data, idx_y+1, idx_x  , gridXsize),
                       LSPI_MATRIX_ELEMENT(data, idx_y  , idx_x+1, gridXsize),
                       LSPI_MATRIX_ELEMENT(data, idx_y+1, idx_x+1, gridXsize)
                      };
    // current_input: [x_cur, y_cur, lower, upper, left, right]
    //           idx:  0      1      2      3      4     5   
    T current_input[6]= {qX, qY, lower, upper, left, right};
    LSPI_VERBOSE_2D(
        displayVec(data_input   , 8, "data_input   ");
        displayVec(current_input, 6, "current_input");
    )
    T result = eval_interpolation_2D(data_input,current_input);
    LSPI_VERBOSE_2D( displayText("Interpolation at p = (%g, %g):  %g\n", qX, qY, result); )
    return result;
} 


} // end of namespace LSPI

#endif // include guard