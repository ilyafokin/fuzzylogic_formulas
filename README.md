The derivation of forumlas for calculating moments of multiplicity distributions using Fuzzy Logic following the paper Phys. Rev. C 110, 064910 by A. Rustamov and a small toy simulation to check the validity.

The derivation is done using sympy in `formulas_derivation.ipynb`. The highest order of the moments and the number of particle species can be set in the second cell. The explicit formulas for the elements of the matrix A and the order of the moments in the vector b can be copied from the output of the cells.

After small adjustments for syntax, the formulas are used in `fuzzy_logic.jl` to calculate moments in a toy simulation for three particles up to fourth order.