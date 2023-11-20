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
  * **Secant**: Even better convergence than the Newton-Raphson Method. In all tests the first iteration managed to estimate a $x$ very close to the solution of the problem. The _mean_ inclination used by the Secant method proved better than the local inclination (derivative) of the Newton-Raphson method.

## Interpolation

* Practical example: **VFP Table Interpolation**
  * Numerical porous-flow models usually rely on VFP tables to associate subsurface and surface flow conditions. Interpolation of the VFP table in the better known simulators is linear.
  * The proposed problem is to compare linear and Spline interpolation with VFP tables.
  * VFP tables from [Unisim-II-H](https://www.unisim.cepetro.unicamp.br/benchmarks/en/unisim-ii/unisim-ii-h) were used as examples. Half of the input variables were ignored when using the interpolation methods and the _missing_ data was compared to the interpolated values.
    * These files were a part of Unisim-II-H when accessed in Nov/19. In the current version (Sep/23) the model no longer uses VFP tables.
  * For each possible combination of parameters in the VFP table, a curve of (LiqFlow, BHP) was extracted. One of the values was ignored, and a spline built with the remaining data. And then the difference between the original BHP and interpolated with the spline was recorded. The same was done with a linear interpolator.
  * Although it is not recommended, some tests used the interpolation methods for extrapolation.
  * The results show that, for the example studied and this particular problem, spline interpolation doesn't offer better results. Linear interpolation is computationally easier and the errors are slightly better than natural spline interpolation.
  * On the second evaluation it was defined a _general rule_ to estimate the first derivative in the extremes of the data points, to be able to use fixed splines. The results improved, and the fixed splines were better than the linear interpolation.

* Methods implemented:
  * Linear interpolation: Most basic for of interpolation (maybe nearest-neighbor can be considered _simpler_). Used as a _benchmark_ to compare with Splines.
  * Splines: Implemented Natural (zero $2^{nd}$ derivative in the extremes) and Fixed Splines (defined $1^{st}$ derivative in the extremes).

## Initial Value Problem

* Practical example: **Fetkovich Aquifer**
  * Compared the computaional efficiency of the Fetkovich Method to solve an analytical aquifer to different methods of solving initial problem values.
  * The first proposed problem was an aquifer connected to a reservoir that was at a different pressure in the initial time.
    * The $4^{th}$ and $5^{th}$ order Runge-Kutta method had results similar to or better than the Fetkovich Method.
    * The multi-steps methods had better performance than Rung-Kutta methods of similar local error.
    * The best results were achieved with the Adams-Bashforth 5 step method (explicit).

* Implementation:
  * Defined a commmon object that holds all methods implemented.
  * Methods implemented:
    * Euler:
      * Can give an error limit (very bad estimate).
      * Option to use Euler together with Aitken to improve responses. The algorithm runs the Euler method three times with increasing number of time-steps ($n$, $2n$, $4n$) and applies the Aitken method in the common times. Accuracy increased by 1-2 orders of magnitude.
    * Runge-Kutta:
      * Implemented $2^{nd}$ (midpoint), $3^{rd}$ (Heu), $4^{th}$ (_original_ version) and $5^{th}$ ($5^{th}$ order Fehlberg) order versions (see [Wikipedia](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods)).
      * Implemented Aitken scheme as well, but results were not so good as with Euler.
    * Multi-steps:
      * Implemented explicit and implicit versions:
        * Adams-Bashforth with 2, 3, 4 and 5 steps (explicit).
        * Adams-Moulton with 1, 2, 3 and 4 steps (implicit).
      * Equivalent Runge-Kutta (same $\mathcal{O}(h^n)$) to find the initial values.
      * Implicit versions use the equivalent explicit version (same $\mathcal{O}(h^n)$) to predict the values (predictor-corrector method).
      * Although it is not the usual method, the possibility to iterate more than once in each time-step was implemented. Default value is just on corrector iteration.
