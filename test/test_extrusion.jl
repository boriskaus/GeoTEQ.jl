using Test, GeoTEQ, ReadVTK


fname_vts_input = "example-data/restraining_bend-cellfields-e2_v2.vts"
fname_vts_output = "extruded.vts"

# extrude the mesh in all directions
extrude_vts(fname_vts_input;
                  fname_vts_output,
                  xmin=3, xmax=5, Δx=1e4,
                  ymin=4, ymax=6, Δy=1e4,
                  zmin=5, zmax=7, Δz=1e4)

# read back data & check size
vtk     = VTKFile(fname_vts_output)
x, y, z     = get_coordinates(vtk)                  
@test size(x) == (265, 75, 141)
