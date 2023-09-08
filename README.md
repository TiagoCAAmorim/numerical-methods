# Numerical Methods
Code for IM253 (Numerical Methods) exercises and final project in C.

Reports (in Latex) in Brazilian Portuguese.

## Root of 1D functions

* Practical example: **Inverse of the Minimum Curvature**
    * Calculation of drilling parameters as a function of cartesian coordinates using de Minimum Curvature Method ([SPE-84246-MS](https://doi.org/10.2118/84246-MS)). This problem is implicit, with the form $g(x)=x$, which can be transformed into $f(x) = g(x) - x = 0$.
    * Developed derivative of the Inverse of the Minimum Curvature for the Newton-Raphson method.
    * Included Mathematica notebook with formulas and additional tests.
    * The best results were found with the Secant model.
* Methods implemented:
  * **Bissection**: All tests converged linearly.
  * **Newton-Raphson**: Convergence was more than quadratic, for $f(x)$ for the proposed problem is closely a straight line. Needed to include limits to $x$ due to the nature of the proposed problem, that doesn't allow for any value of $x$.
  * **Secant**: Even better convergence than the Newton-Raphson Method. In all tests the first iteration managed to estimate an $x$ very close to the solution of the problem. The _mean_ inclination used by the Secant method proved better than the local inclination (derivative) of the Newton-Raphson method.

## Interpolation

* Practical example: **VFP Table Interpolation**
  * Numerical porous-flow models usually relly on VFP tables to associate subsurface and surface flow conditions. Interpolation of the VFP table in the better known simulators is linear.
  * The proposed problem is to compare linear and Spline interpolation with VFP tables.
  * VFP tables from [Unisim-II-H](https://www.unisim.cepetro.unicamp.br/benchmarks/en/unisim-ii/unisim-ii-h) were used as examples. Half of the input variables were ignored when using the interpolation methods and the _missing_ data was compared to the interpolated values.
    * These files were a part of Unisim-II-H when accessed in Nov/19. In the current version (Sep/23) the model no longer uses VFP tables.
  * For each possible combination of parameters in the VFP table, a curve of (LiqFlow,BHP) was extracted. One of the values was ignored and a spline built with the remaing data. And then the difference between the original BHP and interpolated with the spline was recorded. The same was done with a linear interpolator.
  * Although it is not recomended, some of the test used the interpolation methods as extrapolators.
  * The results show that, for the example studied and this particular problem, spline interpolation doesn't offer better results. Linear interpolation is computationally easier and the errors slightly better than spline interpolation.
* Methods implemented:
  * Linear interpolation: Most basic for of interpolation (maybe nearest-neighbour can be considered _simpler_). Used as a _benchmark_ to compare with Splines.
  * Splines: Implemented Natural (zero $2^{nd}$ derivative in the extremes) and Fixed Splines (defined $1^{st}$ derivative in the extremes). For the proposed problem it was used Natural Splines, for there is no information about derivatives.