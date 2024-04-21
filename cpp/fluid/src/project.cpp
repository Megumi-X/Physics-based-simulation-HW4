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
                rhs(idx) = 1.f;
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
    
    ///////////////////////////////////////////////////////////////
}

}
}
