# This file is translated from masbpy.

using LinearAlgebra
using NearestNeighbors
using MeshIO

export MASB, compute_balls!

# FIXME: can't handle duplicate points in input

mutable struct MASB
    D::Dict{String, Vector}
    kdtree::KDTree
    m::Int
    SuperR::Float64
    delta_convergence::Float64
    iteration_limit::Int
    denoise_absmin::Union{Nothing, Float64}
    denoise_delta::Union{Nothing, Float64}
    denoise_min::Union{Nothing, Float64}
    detect_planar::Union{Nothing, Float64}

    function MASB(datadict::Dict{String, Vector}, max_r::Float64; denoise_absmin::Union{Nothing, Float64}=nothing, denoise_delta::Union{Nothing, Float64}=nothing, denoise_min::Union{Nothing, Float64}=nothing, detect_planar::Union{Nothing, Float64}=nothing)
        m = length(datadict["coords"])
        D = Dict(
            "coords" => datadict["coords"],
            "normals" => datadict["normals"],
            "ma_coords_in" => zero(datadict["coords"]).*NaN,
            "ma_coords_out" => zero(datadict["coords"]).*NaN,
            "ma_q_in" => fill(NaN, m),
            "ma_q_out" => fill(NaN, m)
        )
        kdtree = KDTree(datadict["coords"])
        new(D, kdtree, m, max_r, 0.001, 30, denoise_absmin, denoise_delta, denoise_min, detect_planar)
    end
end


function MASB(data::Mesh, max_r::Number; denoise_absmin::Union{Nothing, Float64}=nothing, denoise_delta::Union{Nothing, Float64}=nothing, denoise_min::Union{Nothing, Float64}=nothing, detect_planar::Union{Nothing, Float64}=nothing) 

    datadict = Dict("coords"=> data.position, "normals"=>data.normals)
    return MASB(datadict, Float64(max_r); denoise_absmin=denoise_absmin, denoise_delta=denoise_delta, denoise_min=denoise_min, detect_planar=detect_planar) 

end
    

"""
        compute_balls!(self::MASB)

Compute balls for all points in the dataset.
"""
function compute_balls!(self::MASB)
    for inner in [true, false]
        compute_balls_oneside!(self, inner)
    end
end

"""
     compute_balls_oneside!(self::MASB, inner::Bool=true)
     
Balls shrinking algorithm. Set `inner` to False when outer balls are wanted.
"""
function compute_balls_oneside!(self::MASB, inner::Bool=true)

    # iterate over all point-normal pairs
    for p_i in 1:self.m
        p, n = self.D["coords"][p_i], self.D["normals"][p_i]
        #-- p is the point along whose normal n we are shrinking a ball, its index is p_i
        if !inner
            n = -n
        end

        # initialize some helper variables:
        #-- q will represent the second point that defines a ball together with p and n
        q = nothing
        #-- q_i is the index of q
        q_i = nothing
        #-- r represents the ball radius found in the current iteration (i.e. of the while loop below)
        r = nothing
        #-- r_previous represents the ball radius found in the previous iteration
        r_previous = self.SuperR
        #-- c is the ball's center point in the current iteration
        c = nothing
        #-- c_previous is the ball's center point in the previous iteration
        c_previous = nothing
        #-- j counts the iterations
        j = -1

        while true
            # increment iteration counter
            j += 1
            # set r to last found radius if this isn't the first iteration
            if j > 0
                r_previous = r
            end

            # compute ball center
            c = p .- n .* r_previous

            # keep track of this for denoising purpose
            q_i_previous = q_i

            ### FINDING NEAREST NEIGHBOR OF c

            # find closest point to c and assign to q
            indices, dists  = knn(self.kdtree, c, 2, true)
            candidate_c = self.D["coords"][indices]

            q = candidate_c[1]
            q_i = indices[1]

            # What to do if closest point is p itself?
            if q == p
                # 1) if r_previous==SuperR, apparantly no other points on the halfspace spanned by -n => that's an infinite ball
                if r_previous == self.SuperR
                    r = r_previous
                    break
                else
                    # 2) otherwise just pick the second closest point
                    q = candidate_c[2]
                    q_i = indices[2]
                end
            end

            ### END FINDING NEAREST NEIGHBOR OF c
            # compute new candidate radius r
            r = compute_radius(p, n, q)
            
            ### EXCEPTIONAL CASES

            # if r < 0 closest point was on the wrong side of plane with normal n => start over with SuperRadius on the right side of that plane
            if r < 0
                r = self.SuperR
            # if r > SuperR, stop now because otherwise in case of planar surface point configuration, we end up in an infinite loop
            elseif r > self.SuperR
                r = self.SuperR
                break
            end

            ### END EXCEPTIONAL CASES

            ### DENOISING STUFF
            # i.e. terminate iteration early if certain conditions are satisfied based on (history of) ball metrics

            c_ = p .- n .* r
            # this seems to work well against noisy ma points.
            if !isnothing(self.denoise_absmin)
                if acos(cos_angle(p .- c_, q .- c_)) < self.denoise_absmin && j > 0 && r > norm(q .- p)
                    # keep previous radius:
                    r = r_previous
                    q_i = q_i_previous
                    break
                end
            end

            if !isnothing(self.denoise_delta) && j > 0
                theta_now = acos(cos_angle(p .- c_, q .- c_))
                q_previous = self.D["coords"][q_i_previous]
                theta_prev = acos(cos_angle(p .- c_, q_previous .- c_))

                if (theta_prev - theta_now) > self.denoise_delta && theta_now < self.denoise_min && r > norm(q .- p)
                    # keep previous radius:
                    r = r_previous
                    q_i = q_i_previous
                    break
                end
            end

            if !isnothing(self.detect_planar)
                if acos(cos_angle(q .- p, -n)) > self.detect_planar && j < 2
                    r = self.SuperR
                    break
                end
            end

            ### END DENOISING STUFF

            # stop iteration if r has converged
            if abs(r_previous - r) < self.delta_convergence
                break
            end

            # stop iteration if this looks like an infinite loop:
            if j > self.iteration_limit
                break
            end
        end

        # now store valid points in array (invalid points will be NaN)
        inout = inner ? "in" : "out"

        if r < self.SuperR
            self.D["ma_coords_$inout"][p_i] = c
            self.D["ma_q_$inout"][p_i] = q_i
        end
    end
end