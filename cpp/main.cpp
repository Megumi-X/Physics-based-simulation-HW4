#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"
#include "fluid/include/inviscid_fluid_simulator.hpp"

PYBIND11_MODULE(backend, m) {
    pybind11::class_<backend::fluid::Simulator>(m, "Simulator")
        .def(pybind11::init<const backend::real, backend::Vector2i&, const backend::Vector2i&, const backend::Vector2i&>())
        .def("Forward", &backend::fluid::Simulator::Forward, pybind11::arg("time_step"))
        .def("Project", &backend::fluid::Simulator::Project)
        .def("Advect", &backend::fluid::Simulator::Advect, pybind11::arg("time_step"))
        .def("SolvePoissonEquation", &backend::fluid::Simulator::SolvePoissonEquation, pybind11::arg("div_v"))
        .def("GetVisualizationResult", &backend::fluid::Simulator::GetVisualizationResult, pybind11::arg("resolution"));
}