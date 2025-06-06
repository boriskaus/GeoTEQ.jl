
export compute_median_surface


"""
    xp, Eij_av, Eij_vec = compute_median_surface(ply_fname::String, radius_ma::Float64, radius_cov::Float64)

Computes the median surface from a PLY file containing points and normals.

Input:
- `ply_fname` is the path to the PLY file.
- `radius_ma` is the radius for the moving average smoothing.
- `radius_cov` is the radius for the covariance calculation.

Output:
- `xp` is the array of points on the median surface.
- `Eij_av` is the eigenvector for the minimum eigenvalue, which usually indicates the main orientation of the surface.
- `Eij_vec` is the array of eigenvectors  in x,y,z directions.

Example usage
===
```julia
julia> using 
julia> ply_fname = "example-data/test_surface_1.ply"
julia> radius_ma = 1e4          # Ball radius for moving average    
julia> radius_cov = 25000.0
julia> xp, Eij_av, Eij_vec = compute_median_surface(ply_fname, radius_ma, radius_cov)
```
In a next step, you can filter out some outliers and write the results to a VTK file:
```julia
julia> ind = findall(Eij_av[:,1] .> 0)
julia> write_vts("my_vtp_file", xp[ind], Eij_av[ind,:])
```

"""
function compute_median_surface(ply_fname::String, radius_ma::Number, radius_cov::Number)

    data = load_ply(ply_fname)          # load data
 
    ma = MASB(data, radius_ma)          # Create data structure

    compute_balls!(ma)                  # Perform the computation:

    x   = ma.D["ma_coords_in"]          # Points @ the mid-surface:
    ind = findall(isnan.(x) .== false)  # Extract only valid points
    xp  = x[ind]                        # Valid points  

    Eij_av, Eij_vec, _ = plane_covariance_eigen(xp,radius_cov);

    return xp, Eij_av, Eij_vec
end
