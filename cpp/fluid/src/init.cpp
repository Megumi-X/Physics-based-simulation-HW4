#include "fluid/include/inviscid_fluid_simulator.hpp"

namespace backend {
namespace fluid {

Simulator::Simulator(const real inlet_v, const Vector2i& cell_num, const Vector2i& obstacle_size, const Vector2i& obstacle_bottom_left_corner) :
    cell_num_(cell_num), velocity_x_(cell_num(0) + 1, cell_num(1)), velocity_y_(cell_num(0), cell_num(1) + 1),
    cell_indices_(cell_num(0), cell_num(1)), obstacle_size_(obstacle_size), dof_num_(0), left_inlet_velocity_(inlet_v),
    obstacle_bottom_left_corner_(obstacle_bottom_left_corner) {
    velocity_x_.setZero(); velocity_y_.setZero();
    velocity_x_.block(0, 1, 1, cell_num(1) - 2).fill(inlet_v);
    // Fill in cell_indices_ and dof_num_.
    auto is_inside_obstacle = [&](const integer i, const integer d) -> bool {
        return (i >= obstacle_bottom_left_corner(d)) && (i < obstacle_bottom_left_corner(d) + obstacle_size(d));
    };
    auto invalid_dof = [&](const integer i, const integer j) -> bool {
        return (is_inside_obstacle(i, 0) && is_inside_obstacle(j, 1)) || (j == 0) || (j == cell_num(1) - 1);
    }; 
    for (integer i = 0; i < cell_num(0); ++i)
        for (integer j = 0; j < cell_num(1); ++j)
            if (invalid_dof(i, j)) cell_indices_(i, j) = -1;
            else cell_indices_(i, j) = dof_num_++;
    // Precompute A_inv.
    MatrixXr poisson_matrix(dof_num_, dof_num_);
    // TODO. ////////////////////////////////////////
    poisson_matrix.setZero();
    // For a interior cell, you should fill in something like:
    // poisson_matrix(i, i) = 4.f;
    // poisson_matrix(i, j) = -1.f;
    // (That is, do not multiply any scalars like dx, rho, 1/4 in this matrix.)
    // For boundaries, please refer to the lecture note on how to modify it.
    enum adjacent_type { interior, air, solid };
    for (integer i = 0; i < cell_num(0); ++i) {
        for (integer j = 0; j < cell_num(1); ++j) {
            if (cell_indices_(i, j) == -1) continue;
            for (integer fi = 0; fi < 4; ++fi) {
                integer ni = i + (fi / 2) * ((fi % 2) * 2 - 1);
                integer nj = j + (1 - fi / 2) * ((fi % 2) * 2 - 1);               
                adjacent_type type = interior;
                if (ni >= cell_num(0)) type = air;
                else if (invalid_dof(ni, nj) || ni < 0) type = solid;
                if (type == interior) {
                    poisson_matrix(cell_indices_(i, j), cell_indices_(i, j)) += 1.f;
                    poisson_matrix(cell_indices_(i, j), cell_indices_(ni, nj)) -= 1.f;
                } else if (type == air) {
                    poisson_matrix(cell_indices_(i, j), cell_indices_(i, j)) += 1.f;
                }
            }
        }
    }
    /////////////////////////////////////////////////
    poisson_matrix_inv_ = poisson_matrix.inverse();
}

}
}