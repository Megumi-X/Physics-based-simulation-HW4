#ifndef BACKEND_MAC_GRID_FLUID_SIMULATOR
#define BACKEND_MAC_GRID_FLUID_SIMULATOR

#include "basic/include/config.hpp"
#include "basic/include/log.hpp"

namespace backend {
namespace fluid {

////////////////////////////////////////////////////////////////
// For simplicity, this class only support the following scene.
//
// y
// \  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// |  ->
// |  ->
// |  ->   XXXXXXXX
// |  ->   XX    XX
// |  ->   XX    XX
// |  ->   XXXXXXXX
// |  ->
// |  ->
// |  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// O-----------------------------------------------------> x
//
// There is an inlet flow in the left boundary.
// The top and bottom rows are solid cells.
// There is an axis-aligned square obstacle in the scene.
// The right boundary is a free boundary.
// Each cell has side length 0.1.
// The density is fixed to be 1.0 for simplicity.
// /////////////////////////////////////////////////////////////

class Simulator {
public:
    Simulator(const real inlet_v, const Vector2i& cell_num, const Vector2i& obstacle_size, const Vector2i& obstacle_bottom_left_corner);

    void Forward(const real time_step);
    void Advect(const real time_step);
    void Project();
    const VectorXr SolvePoissonEquation(const VectorXr& div_v) const;
    const std::pair<MatrixXr, MatrixXr> GetVisualizationResult(const integer resolution) const;

private:
    const Vector2r GetVelocity(const Vector2r& position, bool project_to_interior) const;

    const real density_ = 1.f;
    const real dx_ = .1f;
    const real left_inlet_velocity_ = 5.f;
    const Vector2i cell_num_;
    const Vector2i obstacle_size_;
    const Vector2i obstacle_bottom_left_corner_;
    integer dof_num_;
    MatrixXi cell_indices_;
    // Note, v in solid cells should be fixed to be 0.
    MatrixXr velocity_x_;
    // For example, the bottom are solid cells (cells whose y == 0), so
    // both velocity_y_.col(0) (vy stored on the bottom of solid cells)
    // and velocity_y_.col(1) (vy stored on the top of the solid cells) should be 0.
    MatrixXr velocity_y_;
    MatrixXr poisson_matrix_inv_;
};

}
}

#endif