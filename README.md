# GeoTEQ.jl

This is a Julia implementation of the [GeoTEQpy](https://github.com/anthony-jourdon/GeoTEQpy) python package, developed by Anthony Jourdon and coworkers.
It can be used to extract smooth fault planes from geodynamic simulations, where faults form as a result of non-associated plasticity. Typically, such simulations are performed on structured meshes and develop fault zones rather than discrete planes. GeoTEQ can be used, in combination with [ParaView](https://www.paraview.org), to construct smooth fault planes from such data.  


### 1. Installation
Install it in the julia package manager with:
```julia
julia>]
pkg>add https://github.com/boriskaus/GeoTEQ.jl
```
and test it with:
```julia
pkg> test GeoTEQ
``` 

### 2. Differences with the python package
The julia package closely follows the python package, which has a very nice [documentation](https://geoteqpy.readthedocs.io/en/latest/). 
The main differences are:

1. When you extract a fault isocontour from the `xi`, you should save it as a `PLY Polygonal File Format (*.ply)` file, while including info about the normals.
2. The resulting median plane surface is saved as `*.vts` rather than as `*.vtp` file. 

Apart from this, the steps in ParaView are the same. 

### 4. Running the julia package

You can run extract the median plane from the unstructured isosurface with:
```julia
julia> using GeoTEQ
# Specify file names created by saving a surface as *.ply in ParaView:
julia> ply_fname = "example-data/example_contour_strikeslip_1.ply"

# Perform calculations
radius_ma   = 1e4;    # Ball radius
radius_cov  = 25000   # Covariance radius
julia> xp, Eij_av, Eij_vec = compute_median_surface(ply_fname, radius_ma, radius_cov)

# Save output:
julia> write_vts("median_surface", xp, Eij_av)
```


If needed, you can also remove outliers with:
```julia 
julia> ind = findall(Eij_av[:,1].>0)
julia> write_vts("my_vtp_file", xp[ind], Eij_av[ind,:]) # write only active points
```


### 5. Algorithm
The basis of the package is the Shrinking Ball algorithm by Ma et al. (2012), following [this](https://github.com/tudelft3d/masbpy) python implementation, which approximates the Medial Axis Transform (MAT) of an oriented point cloud. We have translated it to julia, for speed reasons and to make it simpler to run on different systems.


