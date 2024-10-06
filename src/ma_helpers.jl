using LinearAlgebra
using NearestNeighbors



"""
    proj(u, v)
"""
function proj(u, v)
    factor = dot(u, v) / dot(u, u)
    return [factor * u[1], factor * u[2], factor * u[3]]
end

"""
    cos_angle(p, q)

Calculate the cosine of angle between vector `p` and `q`
see http://en.wikipedia.org/wiki/Law_of_cosines#Vector_formulation
"""
function cos_angle(p, q)
    cos_angle = dot(p, q) / (norm(p) * norm(q))
    if cos_angle > 1
        return 1
    elseif cos_angle < -1
        return -1
    else
        return cos_angle
    end
end

"""
    compute_radius(p, p_n, q)
compute radius of ball that touches points `p` and `q` and is centered on along the normal `p_n` of `p`.
"""   
function compute_radius(p, p_n, q)
    # this is basic goniometry
    d = norm(p - q)
    cos_theta = dot(p_n, p - q) / d
    return d / (2 * cos_theta)
end


"""
    normalize(vec)
Calculate the normalized vector (norm: one).
"""
normalize(vec) = vec / norm(vec)


function compute_lfs(datadict, k=10)
    # collect all ma_coords that are not NaN
    ma_coords = vcat(datadict["ma_coords_in"], datadict["ma_coords_out"])
    ma_coords = ma_coords[.!any(isnan, ma_coords, dims=2)]

    # build kd-tree of ma_coords to compute lfs
    pykdtree = KDTree(ma_coords)
    if k > 1
        datadict["lfs"] = sqrt.(median(pykdtree.query(datadict["coords"], k)[1], dims=2))
    else
        datadict["lfs"] = sqrt.(pykdtree.query(datadict["coords"], k)[1])
    end
end

"""
    compute_lam(D, inner="in")

Compute for every boundary point p, corresponding ma point m, and other feature point p_ the distance p-p_ 
"""
function compute_lam(D, inner="in")
    D["lam_$inner"] = fill(NaN, size(D["coords"], 1))

    for i in 1:size(D["coords"], 1)
        p = D["coords"][i, :]
        c_p = D["ma_coords_$inner"][i, :]
        if !isnan(c_p[1])
            p_ = D["coords"][D["ma_f2_$inner"][i], :]
            D["lam_$inner"][i] = norm(p - p_)
        end
    end
end

"""
    compute_theta(D, inner="in")
Compute for every boundary point p, corresponding ma point m, and other feature point p_ the angle p-m-p_ 
"""
function compute_theta(D, inner="in")
    D["theta_$inner"] = fill(NaN, size(D["coords"], 1))

    for i in 1:size(D["coords"], 1)
        p = D["coords"][i, :]
        c_p = D["ma_coords_$inner"][i, :]
        if !isnan(c_p[1])
            p_ = D["coords"][D["ma_f2_$inner"][i], :]
            D["theta_$inner"][i] = cos_angle(p - c_p, p_ - c_p)
        end
    end
end
