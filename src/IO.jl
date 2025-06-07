
using WriteVTK, MeshIO, FileIO, ReadVTK
export write_vts, load_ply, extrude_vts

"""
    write_vts(filename::String, xp, Eij_av)

Writes a paraview *.vts file
"""
function write_vts(filename::String, xp, Eij_av)

    vtk_grid(filename, xp) do vtk
        # add datasets...
        vtk["eigv_0"] = Eij_av'
    end
    
    println("Wrote paraview file: $filename.vts")

    return nothing
end


"""
    load_ply(filename::String)
Loads `filename` as a PLY file (generated from a surface in paraview; please save this with surface normals included)
"""
load_ply(filename::String) = load(filename)



"""
    x,y,z, cell_data_names, CellData, point_data_names, PointData = read_vts_file(fname_vtk::String)

Reads a VTS file with cell and point data.
Note that if you have problem reading this data, you may need to open the data in ParaView and save with "save as" 
"""
function read_vts_file(fname_vtk::String)
    vtk     = VTKFile(fname_vtk)
    x, y, z     = get_coordinates(vtk)
    cell_data   = get_cell_data(vtk)
    point_data  = get_point_data(vtk)

    CellData = ();
    cell_data_names = cell_data.names
    for name in cell_data_names
        myCellData = get_data_reshaped(cell_data[name], cell_data = true)
        CellData = (CellData..., myCellData)
    end

    PointData = ();
    point_data_names = point_data.names
    for name in point_data_names
        myPointData = get_data_reshaped(point_data[name], cell_data = false)
        PointData = (PointData..., myPointData)
    end

    return x,y,z, cell_data_names, CellData, point_data_names, PointData
end

"""
    extrude_vts(fname_vts_input::String;
                    fname_vts_output::String="extruded.vts",
                    xmin=nothing, xmax=nothing, Δx=1e4,
                    ymin=nothing, ymax=nothing, Δy=1e4,
                    zmin=nothing, zmax=nothing, Δz=1e4)

Extrudes the vts dataset `fname_vts_input` in `x/y/z` directions by adding rows of data.
The result is saved in `fname_vts_output`.

"""
function extrude_vts(fname_vts_input::String;
                    fname_vts_output="extruded.vts",
                    xmin=nothing, xmax=nothing, Δx=1e4,
                    ymin=nothing, ymax=nothing, Δy=1e4,
                    zmin=nothing, zmax=nothing, Δz=1e4
                    )

    x,y,z, cell_data_names, CellData, point_data_names, PointData = read_vts_file(fname_vts_input)

    for (val, dim, min_row, Δ) in (
        (xmin, 1, true,  Δx),
        (xmax, 1, false, Δx),
        (ymin, 2, true,  Δy),
        (ymax, 2, false, Δy),
        (zmin, 3, true,  Δz),
        (zmax, 3, false, Δz)
    )
        if !isnothing(val)
            x = add_data_rows(x, val, dim; min_row=min_row, Δ=dim==1 ? Δ : 0)
            y = add_data_rows(y, val, dim; min_row=min_row, Δ=dim==2 ? Δ : 0)
            z = add_data_rows(z, val, dim; min_row=min_row, Δ=dim==3 ? Δ : 0)
            CellData  = add_data_rows(CellData,  val, dim; min_row=min_row)
            PointData = add_data_rows(PointData, val, dim; min_row=min_row)
        end
    end

    # save file
    write_vts_file(x,y,z, cell_data_names, CellData, point_data_names, PointData; fname_vtk=fname_vts_output)

    return nothing
end

# adds data rows to the data arrays
function add_data_rows(data::AbstractArray, steps::Int64, dim=2; Δ=0, min_row=true)

    dims = [size(data)...]
    dims[dim] = 1;

    if     dim==1
        data_min = reshape(data[1  ,:,:], dims[1], dims[2], dims[3])        # last column of x
        data_max = reshape(data[end,:,:], dims[1], dims[2], dims[3])        # last column of x
    elseif dim==2
        data_min = reshape(data[:,1  ,:], dims[1], dims[2], dims[3])    # last column of y
        data_max = reshape(data[:,end,:], dims[1], dims[2], dims[3])    # last column of y
    elseif dim==3
        data_min = reshape(data[:,:, 1  ], dims[1], dims[2], dims[3])    # last column of z
        data_max = reshape(data[:,:, end], dims[1], dims[2], dims[3])    # last column of z
    end

    if min_row == true
        for i=1:steps
            data_min .-= Δ
            data = cat(data_min, data; dims=dim)
        end

    elseif min_row == false
        for i=1:steps
            data_max .+= Δ
            data = cat(data, data_max; dims=dim)
        end
    end


    return data
end

# adds data rows for tuples
function add_data_rows(data_tuple::NTuple, steps::Int64, dim=2; min_row=true)

    data_tuple_new = ()
    for data in data_tuple
        data_new = add_data_rows(data, steps, dim; min_row=min_row)
        data_tuple_new = (data_tuple_new..., data_new)
    end

    return data_tuple_new
end


function  write_vts_file(x,y,z, cell_data_names, CellData, point_data_names, PointData; fname_vtk::String="test")

    vtk = vtk_grid(fname_vtk, x, y, z)

    # write cell data
    for i = 1:length(CellData)
        vtk[cell_data_names[i]] = CellData[i]
    end

    # write point data
    for i in 1:length(PointData)
        vtk[point_data_names[i]] = PointData[i]
    end

    close(vtk)
    println("Wrote paraview file: $fname_vtk")
end