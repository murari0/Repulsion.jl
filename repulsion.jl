using LinearAlgebra
using StaticArrays
using Quaternions


function repulsion(npoints, dims, niter = 1e7, convergence_tol::T = 1e-10) where {T <: AbstractFloat}
    R_raw = rand(T, dims, npoints) .- one(T)/2 # generating the first set of random points 
    R = [SVector{dims, T}(ntuple(i -> R_raw[n+i], Val(dims))) for n in 0:dims:(dims*npoints-1)]
    R_n = copy(R) # pre-allocate
    F = zeros(SVector{dims, T},npoints) # pre-allocate "forces"
    for i in 1:niter # starting the conjugate gradient optimization
        F .= zero(T) .* F
        for k in eachindex(R)
            F .-= dist_vector.(R,[R[k]]) # compute the distance in between points
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
	q = Quaternion(p...,true)
	r = rotationmatrix(q)
	# Quaternion to ZYZ Euler angles
	angles[2] = acos(r[3,3])
	if sin(angles[2]) == 0
	    angles[1] = 0
	    angles[3] = asin(r[2,1])
	    if angles[3] < 0
		angles[3] = pi - angles[3]
	    end
	else
	    angles[3] = atan2(r[2,3], r[1,3])
	    if angles[3] < 0
		angles[3] = pi + angles[3]
	    end
	    angles[1] = atan2(r[3,2], -1*r[3,1])
	    if angles[1] < 0
		angles[1] = pi + angles[1]
	    end
	end
    end
    angles
end
