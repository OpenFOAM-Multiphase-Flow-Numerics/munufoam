**[Prerequisites](#Prerequisites)** |
**[Installing](#Installing)** |
**[Documentation](https://openfoam-multiphase-flow-numerics.github.io/numufoam/)** |
# Multiphase


Code repository for the Multiphase Technical Committee


## interface capturing framework

An interface capturing framework with tests and benchmarks


## Getting Started


### Prerequisites

Requires OpenFOAM v2306:

```
https://www.openfoam.com/download/release-history.php
```


### Installing

```
./Allwmake
./get-gmsh.sh # will install gmsh version 4.9.3 as gmshv493
```

### running benchmarks

assuming you have installed python version >= 3.6 (miniconda is a great option)

```
python -m venv env
pip install casefoam

./runAll.sh

or 

python genCases.py
./Allrun
```

### running testsuite (CI)

assuming you have installed python version >= 3.6 (miniconda is a great option)

```
python -m venv env
pip install oftest

py.test # runs the tests
py.test --writeNSteps=1 tutorials/ # test case in run folder
```

### Documentation build

For building the documentation further dependencies like doxygen and sphinx are requirement.

the current  documentation can be fewed [here](https://openfoam-multiphase-flow-numerics.github.io/numufoam/)


## Authors

* **Henning Scheufler**
## License

This project is licensed under the GPL v3 License - see the [LICENSE.md](LICENSE.md) file for details
