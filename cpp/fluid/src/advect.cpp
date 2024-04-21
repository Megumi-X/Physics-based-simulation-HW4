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

    ///////////////////////////////////////////////////
}

}
}