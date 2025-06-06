using MaSB

# Load sample data as ply file (points + normals)
#data = load_ply("example-data/house_dyke_tree.ply")
data = load_ply("example-data/test_surface_1.ply")

# Create data structure
ma = MASB(data, 10)

ma.SuperR = 100 # ball radius

# Perform the computation:
compute_balls!(ma)

# Points @ the mid-surface:
x = ma.D["ma_coords_in"]



# extract only part of them (since we clearly have 2 fault planes)
#x1 = [xl[1] for xl in x];
#ind = findall(x1.>0);

#xp = x[ind]

# 
Eij_vec, Cij_vec = plane_covariance_eigen(xp,1000);

Eij_11 = [Eij[1,1] for Eij in Eij_vec];
Eij_12 = [Eij[1,2] for Eij in Eij_vec];
Eij_13 = [Eij[1,3] for Eij in Eij_vec];
Eij_21 = [Eij[2,1] for Eij in Eij_vec];
Eij_22 = [Eij[2,2] for Eij in Eij_vec];
Eij_23 = [Eij[2,3] for Eij in Eij_vec];
Eij_31 = [Eij[3,1] for Eij in Eij_vec];
Eij_32 = [Eij[3,2] for Eij in Eij_vec];
Eij_33 = [Eij[3,3] for Eij in Eij_vec];




# The medium surface plane is now:
# ma.D["ma_coords_in"]
# 




