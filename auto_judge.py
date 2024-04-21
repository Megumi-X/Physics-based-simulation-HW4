from backend import *
import numpy as np
sim = Simulator(1., [10, 10], [3, 3], [3, 3]) # type: ignore

gt = np.load("test_data.npz")

# Test Poisson.
poisson_ret = sim.SolvePoissonEquation(np.ones(71))
if np.isclose(poisson_ret, gt['p']).all():
    print("Your discretized laplacian operator is consistent with GT.")
else:
    print("Your discretized laplacian operator is possibly wrong.")

# Test projection.
try:
    sim.Project()
    project_ret = sim.GetVisualizationResult(2)
except RuntimeError:
    print("Runtime error occurs in testing projection. The test will exit.")
    exit(0)

if np.isclose(project_ret[0], gt['p0']).all() and np.isclose(project_ret[1], gt['p1']).all():
    print("Your projection is consistent with GT.")
else:
    print("Your projection is possibly wrong.")

# Test advection.
try:
    sim.Advect(0.0001)
    advect_ret = sim.GetVisualizationResult(2)
except RuntimeError:
    print("Runtime error occurs in testing advection. The test will exit.")
    exit(0)

for i in range(advect_ret[0].shape[0]):
    for j in range(advect_ret[0].shape[1]):
        if np.isclose(advect_ret[0][i, j], gt['a0'][i, j]):
            continue
        print(f"Your advection is possibly wrong at position ({i}, {j}).")
        print(f"Yours: {advect_ret[0][i, :]}")
        print(f"GT: {gt['a0'][i, :]}")
for i in range(advect_ret[0].shape[0]):
    for j in range(advect_ret[0].shape[1]):
        if np.isclose(advect_ret[1][i, j], gt['a1'][i, j]):
            continue
        print(f"Your advection is possibly wrong at position ({i}, {j}).")
        print(f"Yours: {advect_ret[1][i, :]}")
        print(f"GT: {gt['a1'][i, :]}")
if np.isclose(advect_ret[0], gt['a0']).all() and np.isclose(advect_ret[1], gt['a1']).all():
    print("Your advection is consistent with GT.")
else:
    print("Your advection is possibly wrong.")
