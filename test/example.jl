using GeoTEQ, Test

# Load sample data as ply file (points + normals)
data = load_ply("example-data/house_dyke_tree.ply")

# Create data structure
ma = MASB(data, 10)

ma.SuperR = 100 # ball radius

# Perform the computation:
compute_balls!(ma)

# Points @ the mid-surface:
x = ma.D["ma_coords_in"]
@test  x[10] ≈ [ -0.542655885219574, 3.6914117336273193, -2.8566672801971436]

# test next dataset
ply_fname = "example-data/test_surface_1.ply"
radius_ma = 100
radius_cov = 50
x1, Eij_av1, Eij_vec1 = compute_median_surface(ply_fname, radius_ma, radius_cov)


@test x1[100] ≈ [  -469.15838623046875, 6563.29931640625, -4168.66650390625]