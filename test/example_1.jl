using GeoTEQ, Test

# load surface with bormal vectors
data = load_ply("example-data/example_contour_strikeslip_1.ply")

@test data[1].points[1][2] ≈ -23701.393f0

# 
radius_ma   = 1e4;    # Ball radius
radius_cov  = 25000   # Covariance radius

ma = MASB(data, radius_ma)          # Create data structure
@test ma.D["coords"][3][2] ≈ -22179.178f0


compute_balls!(ma)                  # Perform the computation:
@test ma.D["ma_q_out"][4] ≈ 80.0

x   = ma.D["ma_coords_in"]          # Points @ the mid-surface:
@test sum(x[1:5]) ≈ [1990.0397033691406, -111891.794921875, 1.065791796875e6]

ind = findall(isnan.(x) .== false)  # Extract only valid points
xp  = x[ind]                        # Valid points  
@test sum(xp) ≈ [   1.8091386696734196e10, -7.745229982509308e8, 9.068905905078125e9]


# Eigenvectors
Eij_av, _ = plane_covariance_eigen(xp,radius_cov);
@test Eij_av[3] ≈ 0.19901120677241035

# filter outliers
ind = findall(Eij_av[:,1].>0)
@test length(ind) == 29947

write_vts("my_vtp_file", xp[ind], Eij_av[ind,:])


# Now do the same with a single function:
ply_fname = "example-data/example_contour_strikeslip_1.ply"
xp1, Eij_av1, Eij_vec1 = compute_median_surface(ply_fname, radius_ma, radius_cov)

# check
@test sum(xp) ≈ sum(xp1)
@test sum(Eij_av) ≈ sum(Eij_av1)