# repulsion.jl
A Julia port of the REPULSION algorithm to generate a uniform grid on a unit hypersphere

## Dependencies
- LinearAlgebra
- StaticArrays
- Quaternions

## Description
This function is a Julia port, by Murari Soundararajan (murari@magnet.fsu.edu), of the REPULSION algorithm originally published by [Bak and Nielsen](https://dx.doi.org/10.1006/jmre.1996.1087) as implemented in Spinach for MATLAB by Ilya Kuprov (https://www.spindynamics.org). Static Arrays are used to improve performance.

    repulsion(npoints, ndims, niter, tol = 1e-10)
    
generates `npoints` number of `ndims` dimensional points and iteratively moves them on the unit hypersphere according to mutual repulsion forces. The algorithm terminates either when `niter` iterations have passed, or when the maximum point displacement between two successive iterations is less than the required tolerance `tol`. The points are returned as a vector of angles that parametrise the points: in two dimensions, a single angle `[0, β, 0]` is returned; in three dimensions, two angles `[0, β, γ]` are returned; and in four dimensions, the points are considered to be quaternions and three ZYZ Euler angles `[α, β, γ]` are returned.

If you use this code, please cite the following papers
- H. J. Hogben, M. Krzystyniak, G. T. P. Charnock, P. J. Hore and I. Kuprov, *Spinach – A software library for simulation of spin dynamics in large spin systems*, J. Magn. Reson., 208 (**2011**) 179–194.
- M. Bak and N. C. Nielsen, *REPULSION, A Novel Approach to Efficient Powder Averaging in Solid-State NMR*, J. Magn. Reson., 125 (**1997**) 132-139. 
