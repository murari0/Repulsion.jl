using LinearAlgebra
using StaticArrays

"""
    repulsion(npoints, dims, niter = 1e7, convergence_tol = 1e-10)

Generates `npoints` random points in `dims` dimensions and distributes them uniformly on the surface of the unit sphere
using the REPULSION algorithm

The function iterates, up to a maximum of `niter` times, until the algorithm converges. The convergence criterion 
`convergence_tol` is interpreted as the target for the maximum displacement between points in two successive iterations.

This function is adapted from its counterpart in Spinach (https://spindynamics.org) for MATLAB. The Spinach 
function, repulsion.m, is Copyright Ilya Kuprov and Frederic Mentink-Vigier
"""
function repulsion(npoints, dims, niter = 1e7, convergence_tol::T = 1e-10) where {T <: AbstractFloat}
    R_raw = rand(T, dims, npoints) .- one(T)/2 # generating the first set of random points 
    R = [SVector{dims, T}(ntuple(i -> R_raw[n+i], Val(dims))) for n in 0:dims:(dims*npoints-1)]
    R_n = copy(R) # pre-allocate
    F = zeros(SVector{dims, T},npoints) # pre-allocate "forces"
    for i in 1:niter # starting the conjugate gradient optimization
        F .= zero(T) .* F
        for k in eachindex(R)
            F .-= dist_vector.(R,[R[k]]) # compute the distance between points
        end
        R_n = normalize.(R_n - dims*F/npoints)
        max_diff = maximum(norm.(R-R_n))
        R = R_n
	    if mod(i, 1000) == 0
	        @debug "Maximum difference " max_diff "At iteration " i # returns the convergence value every 1000 steps
	    end
        if max_diff <= convergence_tol
            @info "Converged at iteration " i # stops and return number of steps when converged
            break
        end
    end
    toangles.(R) # computing the angles
end

"""
    dist_vector(a::SVector{N,T}, b::SVector{N,T}) where N where T

Computes the normalised distance between two vectors `a` and `b`, scaled by the cosine of the angle between them
""" 
function dist_vector(a::SVector{N,T}, b::SVector{N,T}) where N where T
    dd = normalize(a - b)
    if isnan(first(dd))
        zeros(typeof(a))
    else
        dd*dot(a, b)
    end
end

"""
    toangles(p)

Converts a point `p` on the unit sphere to a vector of angles `[α, β, γ]`

The interpretation of `p` and the conversion depend on the number of dimensions:
- Two dimensions: `p` is treated as 2D cartesian co-ordinates, and the polar angle is returned as `β`
- Three dimensions: `p` is treated as 3D cartesian co-ordinates, and the spherical angles `θ` and `φ` are returned as `β` and `γ`
- Four dimensions: `p` is treated as a quaternion, and the ZYZ Euler angles are returned
All angles are shifted to lie in the ranges ``α ∈ [0,2π)``, ``β ∈ [0,π]``, ``γ ∈ [0,2π)``. Missing angles are returned as 0s
"""
function toangles(p::SVector{N,T}) where N where T
    angles = zeros(3)
    if N == 2
	    angles[2] = atan(p[2], p[1])
    elseif N == 3
	    angles[3] = atan(p[2], p[1])
	    angles[2] = pi/2 + atan(p[3], hypot(p[1], p[2]))
    elseif N == 4
	    angles = quat_to_euler(q)
    end
    angles
end

"""
    quat_to_euler(q)

Converts a normalised quaternion to its corresponding Euler angles `[α, β, γ]` about the 
ZYZ rotation axes. 

An active rotation convention is assumed - invert the input quaternion 
for passive rotations. This code is adapted from EasySpin's quat2euler implementation at:
https://github.com/StollLab/EasySpin/blob/main/easyspin/quat2euler.m
EasySpin is Copyright Stefan Stoll and other contributors
"""
function quat_to_euler(q::AbstractVector{T}) where T
    if length(q) != 4
        e = DomainError(q, "A Quaternion must have four elements")
        throw(e)
    end
    
    sy = 2*sqrt((q[3]*q[4] + q[2]*q[1])^2 + (q[2]*q[4] - q[3]*q[1])^2)
    beta = real(atan(sy, 1 - 2*q[2]^2 - 2*q[3]^2))
    if sy < 1e-10 # Singularity check: 1e-10 ≈ 0 for normalised quantities
        alpha = 0.0
        gamma = atan(-2*(q[4]*q[1] - q[2]*q[3]), 1 - 2*(q[2]^2 - q[4]^2))
    else
        alpha = atan(2*(q[3]*q[4] - q[1]*q[2]), 2*(q[2]*q[4] + q[1]*q[3]))
        gamma = atan(2*(q[3]*q[4] + q[1]*q[2]), -2*(q[2]*q[4] - q[1]*q[3]))
        if alpha < 0
            alpha = alpha + 2*pi
        end
        if gamma < 0
            gamma = gamma + 2*pi
        end
    end
    [alpha; beta; gamma]
end