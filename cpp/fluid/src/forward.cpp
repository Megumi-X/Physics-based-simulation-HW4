#include "fluid/include/inviscid_fluid_simulator.hpp"

namespace backend {
namespace fluid {

void Simulator::Forward(const real time_step) {
    Assert(time_step > 0, "fluid::Simulator::Forward", "Non-positive time step.");
    real t = 0;
    bool done = false;
    while (!done) {
        // The 5-dx rule is from Robert Bridson's fluid simulation tutorial.
        real h = 5 * dx_ / (velocity_x_.cwiseAbs().maxCoeff() + velocity_y_.cwiseAbs().maxCoeff());
        if (t + h >= time_step) {
            h = time_step - t;
            done = true;
        }

        // The core steps:
        Advect(h);
        Project();
        t += h;
    }
}

}
}