
We use C++ and Python as the main programming languages. Specifically, we implement the simulation backend in C++ and wrap the interface into a Python package. To achieve this, we require external C++ libraries Eigen 3.4.90 and pybind11 2.11.1, provided in the folder `external`. For visualization, we choose imageio, to be installed by pip. You may refer to the official documentation: https://eigen.tuxfamily.org/ , https://pybind11.readthedocs.io/en/stable/.
## Installation
To run the code, firstly install Python>=3.8. For example in a Conda environment, you could run
```
conda create -n homework4 python=3.10
conda activate homework4
```
Then you should install necessary packages and compile the C++ backend by running `pip install .` in the root directory of our codebase (where there is a `setup.py`). Whenever you change any C++ code, remember to rerun `pip install .` before you run any Python scripts.

## Trouble-shooting
### On Windows or MacOS
We do not officially guarantee that the installation on a Windows or MacOS system is bugfree for every homework, so prepare a Ubuntu environment is recommended for this course. For this homework, if you want a video output, you may find your own way to generate a video from images depending on your platform.

### Type annotations
If you would like to automatically generate the type annotation for the backend package, you may try the following commands in the root directory of the code base:
```
pip install pybind11-stubgen
pybind11-stubgen backend -o=. --numpy-array-remove-parameters --ignore-all-errors
```
It will automatically generate a file `backend.pyi` that explains the exposed APIs to your IDLE.
