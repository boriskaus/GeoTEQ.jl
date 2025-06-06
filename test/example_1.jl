using MaSB

# check how to load the contour with 
data = load_ply("example-data/example_contour_strikeslip_1.ply")


radius_ma = 1e4;    # Ball radius
radius_cov = 25000  # Covariance radius


ma = MASB(data, radius_ma)          # Create data structure
compute_balls!(ma)                  # Perform the computation:
x   = ma.D["ma_coords_in"]          # Points @ the mid-surface:
ind = findall(isnan.(x) .== false)  # Extract only valid points
xp  = x[ind]            # Valid points  

# 
Eij_av, _ = plane_covariance_eigen(xp,radius_cov);


# filter outliers
ind = findall(Eij_av[:,1].>0)

write_vts("my_vtp_file", xp[ind], Eij_av[ind,:])
