# DynamicOptimalTransport

Source code for paper [A Novel Skip Orthogonal List for Dynamic Optimal Transport Problem](https://xyxu2033.github.io/pub/2023_dynamic_ot.pdf)

## Install Dependencies

In order to build our program and run our code, you should first install CMake (>= 3.12) and Python 3 (>= 3.6). It is recommended to install pip.

```bash
sudo apt update && sudo apt install -y cmake python3 python3-pip
```

We use pybind11 library to build our code, and we use functions in Python Optimal Library (POT) as a baseline to measure our performance. It is recommended to run the following line in a virtual environment such as [Conda](https://www.anaconda.com/).

```bash
python3 -m pip install "pybind11[global]" pot
```

## Build and Install Our Library

To build our project through CMake, simply switch to the root directory of our project and execute the following

```bash
mkdir -p build && cd build && cmake .. && make
```

To install it, continue executing the following

```bash
make install
```

## Core Functions

* `tree = make_emd(points, flows, dist_func)` builds a new dynamic optimal transport system. Points have to be `np.array` instances in shape `(n, d)` respectively, where `n` denotes the number of points and `d` denotes the dimension of the space. Flows have to be a 1D `np.array` of length `n`, where the sum of all entries have to be exactly `0`. If `flows[i] > 0`, then point `i` is a supply point, otherwise `i` is a demand point. `dist_func` should be a callable object that takes two 1D `np.array` with length `d` and return their distance. All types have to be integers, i.e., `np.int32` or `np.int64`. The dynamic optimal transport system is stored in `tree`. Points are labeled from `0` to `n - 1` respectively.
* `tree.modify(point_label, new_position)` moves the point labeled with `point_label` to a new point `new_position` in the current space $\mathbb R^d$.
* `tree.sendflow(supply_label, demand_label, flow)` sends `flow` amount of flow from point labeled as `supply_label` to point labeled as `demand_label`.
* `new_label = tree.insert(new_position)` inserts a new point into the system, and it is labeled as `new_label`.
* `tree.delete(label)` deletes the point labeled as `label` from the current system.
