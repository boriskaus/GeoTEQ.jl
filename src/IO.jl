
using WriteVTK
export write_vts, load_ply

"""
    write_vts(filename::String, xp, Eij_av)

Writes a paraview file
"""
function write_vts(filename::String, xp, Eij_av)

    vtk_grid(filename, xp) do vtk
        # add datasets...
        vtk["Eij_av"] = Eij_av'
    end
    
    println("Wrote paraview file: $filename.vts")

    return nothing
end


"""
    load_ply(filename::String)
Loads `filename` as a PLY file (generated from a surface in paraview; please save this with surface normals included)
"""
load_ply(filename::String) = load(filename)
