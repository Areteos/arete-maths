# arete-maths
A collection of mathematical and statistical functions and objects that I've needed to write over the years. Requires the repository [Areteos/arete-utilities](https://github.com/Areteos/arete-utilities) to function, and provides more in-depth functionality than the MathsFunctions class in that repository.

Currently consists of the calculus, linearAlgebra, and kernelDensityEstimation modules.

Unit testing is not yet finished for the entire project, stay tuned!

Every complete item has JavaDoc explaining its use.

## calculus
### differentiation
A package full of effectively immutable Types representing infinitely differentiable elementary mathematical expressions like constants, summation, multiplication, exponentiation, and so on. These can in turn be composed together in a variety of ways to yield more complex expressions: for instance, 

5(3-x)<sup>2</sup>

can be expressed as 5x<sup>2</sup> composed with (3-x).

These expressions are constructed through factory functions, which automatically expand and simplify in an attempt to ensure that any two mathematically equivalent expressions have the same internal representation. Likewise, any mathematical operation on these expressions always attempts to return the simplest or most expanded possible result:

5x<sup>2</sup> times 3x<sup>4</sup> will become 15x<sup>6</sup> rather than 5x<sup>2</sup>3x<sup>4</sup>

This together helps prevent long chains of operations and compositions from resulting in exponentially growing expressions. There will never be a situation for instance where repeat addition results in something like:

(1 + (1 + (1 + (1 + (1 + (1 + (1 + (1 + (1 + x))))))))), this will always flatten and simplify to (9 + x)

These expressions, in addition to allowing for convenient algebraic operations, also allow for algebraic differentiation without resorting to expensive and/or inaccurate numerical approximations. Differentiating a DifferentiableFunction (the parent class of all the expressions described above) yields another DifferentiableFunction. This is done laziy, using a hidden static map of each previously differentiated function to its derivative, which can only be accessed (read or write) via private synchronised static methods.

All this makes for a fast, accurate, and convenient library for algebra and differentiation. 

### differentialEquations
This package is for the purpose of solving partisal differential equations using numerical methods. Currently only contains classes for solving using the finite difference method: an abstract parent class for finite difference solvers in general, and a variably-implicit concrete implmentation. Currently solves up to second-order PDEs with known time- and space-variable coefficients, although work needs to be done to allow for the full range of second-order PDEs.

## linearAlgebra
Contains a tridiagonal and guassian algorithm for solving systems of linear equations. Both can handle systems of arbitrary size, but cannot solve systems that are under-constrained.

## kernelDensityEstimation
### Kernels
Small helper class containing static methods for generating kernels of various bandwidth and shape. Currently implemented using DifferentiableFunctions, since the original focus was on the (infinitely differentiable) Gaussian kernel, but in the future will contain methods for generating undifferentiable kernels.

### KernelDensityEstimator
Class that accepts a list of samples and weights, as well as kernel function, and creates an estimate of the probability density function that generated the original points. As you might expect, uses [kernel density estimation](https://en.wikipedia.org/wiki/Kernel_density_estimation).

### KDEViaDiffusion
Experimental! Uses the techniques described in [this paper](https://www.researchgate.net/publication/47744359_Kernel_Density_Estimation_via_Diffusion) by Zdravko et al to create more advanced kernel density estimations. The `getGaussianEstimator` function is working as intended, and produces a gaussian kernel density estimator with an optimal bandwidth for the given data set. Unfortunately the `getDiffusionEstimator` function is not yet working properly, and I can't recommend using it yet!
