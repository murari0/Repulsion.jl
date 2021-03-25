using LinearAlgebra
using StaticArrays

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

function dist_vector(a::SVector{N,T}, b::SVector{N,T}) where N where T
    dd = normalize(a - b)
    if isnan(first(dd))
        zeros(typeof(a))
    else
        dd*dot(a, b)
    end
end

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