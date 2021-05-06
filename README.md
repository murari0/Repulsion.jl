# Repulsion.jl
A Julia port of the REPULSION algorithm to generate a uniform grid on a unit hypersphere

## Dependencies
- LinearAlgebra
- StaticArrays
- Quaternions

## Optional Dependencies - writing data to binary files
- [MAT.jl](https://github.com/JuliaIO/MAT.jl)
  - To write data as a MATLAB `.mat` file
- [JLD.jl](https://github.com/JuliaIO/JLD.jl)
  - To write data in the Julia Data format `.jld`, an implementation of [HDF5](http://www.hdfgroup.org/HDF5/)
- [JLD2.jl](https://github.com/JuliaIO/JLD2.jl)
  - A native Julia implementation of a subset of HDF5
    
## Description
This function is a Julia port, by Murari Soundararajan (murari@magnet.fsu.edu), of the REPULSION algorithm originally published by [Bak and Nielsen](https://dx.doi.org/10.1006/jmre.1996.1087) as implemented in Spinach for MATLAB by Ilya Kuprov and Frederic Mentink-Vigier (https://www.spindynamics.org). Static Arrays are used to improve performance.

    repulsion(npoints, ndims, niter = 1e7, tol = 1e-10, fileopts...)
    
generates `npoints` number of `ndims` dimensional points and iteratively moves them on the unit hypersphere according to mutual repulsion forces. The algorithm terminates either when `niter` iterations have passed, or when the maximum point displacement between two successive iterations is less than the required tolerance `tol`. The points are returned as a vector of angles that parametrise the points: in two dimensions, a single angle `[0, β, 0]` is returned; in three dimensions, two angles `[0, β, γ]` are returned; and in four dimensions, the points are considered to be quaternions and three ZYZ Euler angles `[α, β, γ]` are returned.

The resulting grid is written to a file specified by the optional `fileopts` arguments. If present, these must be in the form of up to two additional String arguments. The first is a case-insensitive file format specifier, and the second is an optional filename string. If the filename is absent, a default file name is generated from the current system time. The following file formats are currently supported:
- "CSV": the grid is written as a comma separated list of angles, with each point on a new line
- "TSV": the grid is written as a tab separated list of angles, with each point on a new line
- "MAT": (Requires MAT.jl) the grid is separated into four individual arrays of length `npoints` holding the alpha, beta, gamma angles and weights associated with each point. The four arrays, named `alphas`, `betas`, `gammas` and `weights` are saved in the file
- "JLD"/"JLD2": (Requires JLD.jl / JLD2.jl) the grid is saved as the variable `angles`, a vector of `npoints` where each element is a 4-element vector holding the alpha, beta, gamma angles and weights of each point

If an unsupported file format is specified, or if the file format requires an external package that has not already been loaded by the user, the data is written to TSV as a fallback

## Citing
If you use this code, please cite the following papers
- H. J. Hogben, M. Krzystyniak, G. T. P. Charnock, P. J. Hore and I. Kuprov, *Spinach – A software library for simulation of spin dynamics in large spin systems*, J. Magn. Reson., 208 (**2011**) 179–194.
- M. Bak and N. C. Nielsen, *REPULSION, A Novel Approach to Efficient Powder Averaging in Solid-State NMR*, J. Magn. Reson., 125 (**1997**) 132-139. 
