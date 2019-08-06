# Vortex_in_Cell_Method
Implementation of a Two Dimensional Vortex-in-Cell Method for Fluid Dynamics: This code documents the process I went through to gain an understanding of fluid dynamics using a two-dimensional vortex particle in cell method. Utilizing ordinary and partial differential equations, we took an in-depth look at the interaction between particles based off their vorticity strength in two dimensions. This code later examines three tests to determine the effectiveness.

## Getting Started

This code has the option to use two different solvers to solve for Poissons equation. File ConvSolver.py is the default that we will be using, if you would like to use the other file PoissonSolver which is a C++ wrapper follow this [link]() for installing the required packages. The following steps will assume you are working with the ConvSolver.py. 

### Prerequisites

This code assumes that you have and can run Python3, if you do not follow this [link](https://realpython.com/installing-python/)

## Running the tests

All the test should run on their own. In order to avoid confusion when reading the vorticity strength of each particle they have been left out of the test sections below and are better explained in the supporting PDF files.

### Test 1

Test one included one particle placed near the center of the grid. Any self-forcing is nonphysical so over time we should not expect to see this particle move.

### Test 2

Test two contained two particles, one particle at (0.25, 0.5) and another at (0.75, 0.5). The zero-strength particle should rotate around the other particle while the other particle in the center stays fixed.

### Test 3

Test three involved two particles and allowed us to see if they interact correctly. One particle is located at (0.25, 0.5) and another at (0.75, 0.5). The expectation is that both particles should rotate around the center (0.5, 0.5).

### Test 4

Test four is best explained looking at the supporting PDF document titled ""

## Deployment

Assuming you have changed into the directory that contains the both the Vortex_in_Cell.py file and the ConvSolver.py file, in your terminal window run the following code:

```
python3 Vortex_in_Cell.py
```

## Authors

* **David Akinyemi** - Bates College - Summer 2019

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thank you to Colin Wahl for creating both the PoissonSolver and ConvSolver.
