# GeoTEQ.jl
[![CI](https://github.com/boriskaus/GeoTEQ.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/boriskaus/GeoTEQ.jl/actions/workflows/CI.yml)
[![DOI](https://zenodo.org/badge/868229553.svg)](https://doi.org/10.5281/zenodo.15645783)

This is a Julia implementation of the [GeoTEQpy](https://github.com/anthony-jourdon/GeoTEQpy) python package, developed by Anthony Jourdon and coworkers.
It can be used to extract smooth fault planes from geodynamic simulations, where faults form as a result of non-associated plasticity. Typically, such simulations are performed on structured meshes and develop fault zones rather than discrete planes. GeoTEQ can be used, in combination with [ParaView](https://www.paraview.org), to construct smooth fault planes from such data.  

![strike_slip_example](./doc/img/strike_slip.png)
![bend_rift_example](./doc/img/bend_rift.png)

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

Or extrude the original `*.vts` file:
```julia 
julia> fname_vts_input = "example-datarestraining_bend-cellfields-e2_v2.vts"

# extrude the mesh in y-direction:
julia> extrude_vts(fname_vts_input; fname_vts_output="extruded.vts", ymin=4, ymax=6, Δy=1e4)
```


### 5. Algorithm
The basis of the package is the Shrinking Ball algorithm by Ma et al. (2012), following [this](https://github.com/tudelft3d/masbpy) python implementation, which approximates the Medial Axis Transform (MAT) of an oriented point cloud. We have translated it to julia, for speed reasons and to make it simpler to run on different systems.


### 6. Funding
Funding to develop this julia interface was provided by the [ChEESE-2p](https://cheese2.eu) project. 
