# 1. Installation and Uninstallation

## 1.1 Installation

General Installation Method:

```bash
python3 scripts/mk_make.py
cd build
make
sudo make install
```

You can select the installation path:

```bash
python3 scripts/mk_make.py --prefix=/home/leo
cd build
make
make install
```

Installation can be done in `debug` mode.

```bash
python3 scripts/mk_make.py --debug --trace
cd build
make
sudo make install
```

See `mk_make.py` for optional commands:

`python3 scripts/mk_make.py -h`

## 1.2 Uninstallation

```bash
cd build
sudo make uninstall
```

Cleans up an environment that has been polluted by python scripts:

`git clean -fx src`

## 1.3 Using Z3 in a virtual environment

Referencing the z3 package in python

```shell
python3 scripts/mk_make.py --python
cd build
make
make install
```

# 2. Optional trace options

```bash
Z3_enable_trace("assumptions");
```

# 3. Adding Components

Add Cmake file CMakeList.txt.

e.g.

```cmake
z3_add_component(sampler
  SOURCES
    sampler.cpp
)
```

If there is a dependency between the components to be added, the dependency needs to be defined in `scripts/mk_project.py`.

e.g.

```python
add_lib('sampler', ['util'], 'sampler')
```

The function used is:

```python
def add_lib(name, deps=[], path=None, includes2install=[]):
    c = LibComponent(name, path, deps, includes2install)
    reg_component(name, c)
```

note: Remember to use a python script to recompile the Makefile after modifying the CMake file.