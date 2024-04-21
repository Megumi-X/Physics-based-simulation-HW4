#include "fluid/include/inviscid_fluid_simulator.hpp"

namespace backend {
namespace fluid {

void Simulator::Advect(const real time_step) {
    auto RK2 = [&] (const Vector2r& position, const real h) -> const Vector2r  {
        Vector2r velocity_now = GetVelocity(position, true);
        Vector2r velocity_half = GetVelocity(position + 0.5f * h * velocity_now, true);
        return position + h * velocity_half;
    };
    // TODO. Apply semi-Lagrange advection with RK-2.//
    MatrixXr new_velocity_x = velocity_x_;
    MatrixXr new_velocity_y = velocity_y_;
    for (integer i = 0; i < cell_num_(0) + 1; ++i) {
        for (integer j = 0; j < cell_num_(1); ++j) {
            if (i == 0) continue;
            if ((cell_indices_(i, j) == -1 || cell_indices_(i - 1, j) == -1) && i != cell_num_(0)) continue;
            Vector2r position = Vector2r{static_cast<real>(i), static_cast<real>(j) + 0.5f} * dx_;
            Vector2r new_position = RK2(position, -time_step);
            new_velocity_x(i, j) = GetVelocity(new_position, true)(0);
        }
    }
    for (integer i = 0; i < cell_num_(0); ++i) {
        for (integer j = 0; j < cell_num_(1) + 1; ++j) {
            if (j == 0 || j == cell_num_(1)) continue;
            if (cell_indices_(i, j) == -1 || cell_indices_(i, j - 1) == -1) continue;
            Vector2r position = Vector2r{static_cast<real>(i) + 0.5f, static_cast<real>(j)} * dx_;
            Vector2r new_position = RK2(position, -time_step);
            new_velocity_y(i, j) = GetVelocity(new_position, true)(1);
        }
    }
    velocity_x_ = new_velocity_x;
    velocity_y_ = new_velocity_y;
    ///////////////////////////////////////////////////
}

}
}