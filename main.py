import backend
import numpy as np
import imageio

nx = 60
ny = 20
resolution = 20
inlet_v = 5.
sim = backend.Simulator(inlet_v, [nx, ny], [6, 6], [7, 7])
sim.Project()
for i in range(600):
    mat_x, mat_y = sim.GetVisualizationResult(resolution)
    max_x = np.max(mat_x)
    mat_arr = np.zeros((nx * resolution, ny * resolution, 3), dtype=np.float32)
    mat_arr[:, :, 1] = mat_x / inlet_v / 4 + 0.5
    mat_arr[:, :, 0] = 0.1
    mat_arr[:, :, 2] = mat_y / inlet_v / 4 + 0.5
    mat_arr[(mat_x ** 2 + mat_y ** 2) < 1e-10] = 0.
    mat_arr *= 255
    imageio.imwrite(f"output/{str(i).zfill(3)}.png", np.swapaxes(mat_arr.astype(np.uint8), 0, 1))
    sim.Forward(1 / 30)

# You may also generate the video if your system supports it.
# writer = imageio.get_writer('output.mp4', fps=30)
# for file_name in [f"output/{str(i).zfill(3)}.png" for i in range(600)]:
#     writer.append_data(imageio.imread(file_name))
# writer.close()