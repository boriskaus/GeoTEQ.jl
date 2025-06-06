# This routine follows the idea of Jourdon et al. (2024) to smoothen the fault plane:
# https://arxiv.org/pdf/2407.20609
# 
using NearestNeighbors, Statistics, LinearAlgebra

export plane_covariance_eigen

"""
    Eij_av, Eij_vec, Eij_val, Cij = plane_covariance_eigen(x::Vector{Point{N,_T}}, radius) where {N,_T}

Computes the covariance matrix `Cij`, eigenvectors `Eij_vec`, and eigenvalues `Eij_val` of points of a medium axes plane that are within a given radius.
`Eij_av` is the eigenvector corresponding to the smallest eigenvalue, which represents the average direction of the medium surface.
"""
function plane_covariance_eigen(x::Vector{Point{N,_T}}, radius) where {N,_T}

    kdtree = BallTree(x)

    Cij  = zeros(_T, N,N)
    Cij_vec     = fill(Cij,length(x))
    Eij_vec     = fill(Cij,length(x))
    Eij_vec_x   = zeros(length(x),3)
    Eij_vec_y   = zeros(length(x),3)
    Eij_vec_z   = zeros(length(x),3)
    Eij_av      = zeros(length(x),3)
    
    Eij_val_x   = fill(0.0,length(x))
    Eij_val_y   = fill(0.0,length(x))
    Eij_val_z   = fill(0.0,length(x))

    for id in eachindex(x) 
        # find all points that are within the radius of the point
        idxs  = inrange(kdtree, x[id], radius)
        x_pts = x[idxs];
        x_m   = mean(x_pts)
      
        npts = length(x_pts)
        Cij  = zeros(_T, N,N)
        for i=1:N
            for j=1:N
                for pt = 1:npts
                    Cij[i,j]  += (x_pts[pt][i] - x_m[i])*(x_pts[pt][j] - x_m[j])
                end
                Cij[i,j] /= npts
            end
        end
       
        Cij_vec[id] = Cij       # covariance matrix of each point
        Eij = -eigvecs(Cij)      # eigenvectors
        A   = eigvals(Cij)      # eigenvalues
        
       

        Eij_vec[id] = Eij
        Eij_vec_x[id,:] = Eij[:,1]
        Eij_vec_y[id,:] = Eij[:,2]
        Eij_vec_z[id,:] = Eij[:,3]

        _,ind = findmin(A)         # smallest eigenvalue index
        Eij_av[id,:] = Eij[:,ind] # eigenvector corresponding to the smallest eigenvalue

        Eij_val_x[id] = A[1]
        Eij_val_y[id] = A[2]
        Eij_val_z[id] = A[3]
    end



    return Eij_av, (Eij_vec_x, Eij_vec_y, Eij_vec_z), (Eij_val_x, Eij_val_y, Eij_val_z), Cij_vec, Eij_vec
end