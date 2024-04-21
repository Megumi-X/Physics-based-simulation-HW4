#include "fluid/include/inviscid_fluid_simulator.hpp"

namespace backend {
namespace fluid {

const VectorXr Simulator::SolvePoissonEquation(const VectorXr& div_v) const {
    return poisson_matrix_inv_ * div_v;
}

void Simulator::Project() {
    VectorXr rhs = VectorXr::Zero(dof_num_);
    for (integer i = 0; i < cell_num_(0); ++i)
        for (integer j = 0; j < cell_num_(1); ++j) {
            const integer idx = cell_indices_(i, j);
            if (idx != -1) {
                // TODO. What is the right-hand-size of the Poisson's equation?
                real v_x_0 = velocity_x_(i, j);
                real v_x_1 = velocity_x_(i + 1, j);
                real v_y_0 = velocity_y_(i, j);
                real v_y_1 = velocity_y_(i, j + 1);
                rhs(idx) = v_x_1 - v_x_0 + v_y_1 - v_y_0;
                ///////////////////////////////////////////////////////////////
            }
        }
    const VectorXr sol = SolvePoissonEquation(rhs);
    MatrixXr pressure_omit_scalars(cell_num_(0), cell_num_(1));
    pressure_omit_scalars.setZero();
    for (integer i = 0; i < cell_num_(0); ++i)
        for (integer j = 0; j < cell_num_(1); ++j) {
            const integer idx = cell_indices_(i, j);
            if (idx != -1) pressure_omit_scalars(i, j) = sol(idx);
        }

    // TODO. Update the velocity according to the pressure.
    for (integer i = 0; i < cell_num_(0) + 1; ++i) {
        for (integer j = 0; j < cell_num_(1); ++j) {
            if (i == 0) continue;
            if ((cell_indices_(i, j) == -1 || cell_indices_(i - 1, j) == -1) && i != cell_num_(0)) continue;
            if (i < cell_num_(0))
                velocity_x_(i, j) += sol(cell_indices_(i, j)) - sol(cell_indices_(i - 1, j));
            else if (i == cell_num_(0))
                velocity_x_(i, j) += -sol(cell_indices_(i - 1, j));
        }
    }
    for (integer i = 0; i < cell_num_(0); ++i) {
        for (integer j = 0; j < cell_num_(1) + 1; ++j) {
            if (j == 0 || j == cell_num_(1)) continue;
            if (cell_indices_(i, j) == -1 || cell_indices_(i, j - 1) == -1) continue;
            velocity_y_(i, j) += sol(cell_indices_(i, j)) - sol(cell_indices_(i, j - 1));
        }
    }
    ///////////////////////////////////////////////////////////////
}

}
}
