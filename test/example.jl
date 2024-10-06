using FileIO
using MaSB

# Load sample data as ply file (points + normals)
#data = load("example-data/house_dyke_tree.ply")
data = load("example-data/test_surface_1.ply")

# Create data structure
ma = MASB(data, 10)

ma.SuperR = 100 # ball radius

# Perform the computation:
compute_balls!(ma)


# The medium surface plane is now:
# ma.D["ma_coords_in"]
# 




