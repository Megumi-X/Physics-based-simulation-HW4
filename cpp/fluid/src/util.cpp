#include "fluid/include/inviscid_fluid_simulator.hpp"

namespace backend {
namespace fluid {

static const real Interpolate(const real lo, const real hi, const real w) { return w * hi + (1 - w) * lo; }

const Vector2r Simulator::GetVelocity(const Vector2r& position, bool project_to_interior) const {
    const real eps = 1e-7;
    if (position(0) < eps) return {left_inlet_velocity_, 0.0};
    Vector2r lattice_position = position / dx_;
    if (project_to_interior) {
        if (lattice_position(0) > cell_num_(0) - eps)
            lattice_position(0) = cell_num_(0) - eps;
        if (lattice_position(1) < 1 + eps)
            lattice_position(1) = 1 + eps;
        if (lattice_position(1) > cell_num_(1) - 1 - eps)
            lattice_position(1) = cell_num_(1) - 1 - eps;

        // Check if within obstacle bounds
        Eigen::Index idx;
        if (Vector4r{
            lattice_position(0) - obstacle_bottom_left_corner_(0),
            obstacle_bottom_left_corner_(0) + obstacle_size_(0) - lattice_position(0),
            lattice_position(1) - obstacle_bottom_left_corner_(1),
            obstacle_bottom_left_corner_(1) + obstacle_size_(1) - lattice_position(1)
        }.minCoeff(&idx) > -eps) std::array<std::function<void(void)>, 4>{
            [&] { lattice_position(0) = obstacle_bottom_left_corner_(0) - eps; },
            [&] { lattice_position(0) = obstacle_bottom_left_corner_(0) + obstacle_size_(0) + eps; },
            [&] { lattice_position(1) = obstacle_bottom_left_corner_(1) - eps; },
            [&] { lattice_position(1) = obstacle_bottom_left_corner_(1) + obstacle_size_(1) + eps; }
        }[idx]();
    }
    auto vx = [&] {
        const integer i = static_cast<integer>(lattice_position(0));
        const integer j = static_cast<integer>(lattice_position(1) - 0.5f);
        Assert(i >= 0 && i <= cell_num_(0), "backend::fluid::Simulator::GetVelocity::vx", "Invalid inputs. You may need to use project_to_interior=true.");
        const Vector2r weights = lattice_position - Vector2r{ static_cast<real>(i), static_cast<real>(j) + 0.5f };
        if (j == 0) return Interpolate(velocity_x_(i, 1), velocity_x_(i + 1, 1), weights(0));
        if (j == cell_num_(1) - 2) return Interpolate(velocity_x_(i, j), velocity_x_(i + 1, j), weights(0));
        if (j == obstacle_bottom_left_corner_(1) - 1 && i >= obstacle_bottom_left_corner_(0) && i < obstacle_bottom_left_corner_(0) + obstacle_size_(0))
            return Interpolate(velocity_x_(i, j), velocity_x_(i + 1, j), weights(0));
        if (j == obstacle_bottom_left_corner_(1) + obstacle_size_(1) - 1 && i >= obstacle_bottom_left_corner_(0) && i < obstacle_bottom_left_corner_(0) + obstacle_size_(0))
            return Interpolate(velocity_x_(i, j + 1), velocity_x_(i + 1, j + 1), weights(0));
        const real x0 = Interpolate(velocity_x_(i, j), velocity_x_(i, j + 1), weights(1));
        const real x1 = Interpolate(velocity_x_(i + 1, j), velocity_x_(i + 1, j + 1), weights(1));
        return Interpolate(x0, x1, weights(0));
    };
    auto vy = [&] {
        const integer i = static_cast<integer>(lattice_position(0) - 0.5f);
        const integer j = static_cast<integer>(lattice_position(1));
        Assert(j >= 1 && j < cell_num_(1), "backend::fluid::Simulator::GetVelocity::vy", "Invalid inputs. You may need to use project_to_interior=true.");
        const Vector2r weights = lattice_position - Vector2r{ static_cast<real>(i) + 0.5f, static_cast<real>(j) };
        if (i == -1) return Interpolate(0.0, Interpolate(velocity_y_(0, j), velocity_y_(0, j + 1), weights(1)), weights(0));
        if (i == cell_num_(0) - 1) return Interpolate(velocity_y_(i, j), velocity_y_(i, j + 1), weights(1));
        if (i == obstacle_bottom_left_corner_(0) - 1 && j >= obstacle_bottom_left_corner_(1) && j < obstacle_bottom_left_corner_(1) + obstacle_size_(1))
            return Interpolate(velocity_y_(i, j), velocity_y_(i, j + 1), weights(1));
        if (i == obstacle_bottom_left_corner_(0) + obstacle_size_(0) - 1 && j >= obstacle_bottom_left_corner_(1) && j < obstacle_bottom_left_corner_(1) + obstacle_size_(1))
            return Interpolate(velocity_y_(i + 1, j), velocity_y_(i + 1, j + 1), weights(1));
        const real y0 = Interpolate(velocity_y_(i, j), velocity_y_(i, j + 1), weights(1));
        const real y1 = Interpolate(velocity_y_(i + 1, j), velocity_y_(i + 1, j + 1), weights(1));
        return Interpolate(y0, y1, weights(0));
    };
    return { vx(), vy() };
}

const std::pair<MatrixXr, MatrixXr> Simulator::GetVisualizationResult(const integer resolution) const {
    MatrixXr convas_x(cell_num_(0) * resolution, cell_num_(1) * resolution);
    convas_x.setZero();
    MatrixXr convas_y(convas_x);
    for (integer i = 0; i < cell_num_(0); ++i)
        for (integer j = 0; j < cell_num_(1); ++j)
            if (cell_indices_(i, j) != -1)
                for (integer k = 0; k < resolution; ++k)
                    for (integer l = 0; l < resolution; ++l) {
                        const Vector2r v = GetVelocity(Vector2r{ i * dx_ + k * dx_ / resolution, j * dx_ + l * dx_ / resolution }, false);
                        convas_x(i * resolution + k, j * resolution + l) = v(0);
                        convas_y(i * resolution + k, j * resolution + l) = v(1);
                    }
    return { convas_x, convas_y };
}

}
}