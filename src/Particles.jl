
module Particles

using NetCDF

include("unstructured_grid.jl") #support for fast indexing of unstructured grids

include("grids.jl") #generic grid structure

include("cartesian_grid.jl") #support for fast indexing of cartesian (2D) grids

include("wms_client.jl") #support for background images downloaded with WMS

include("dflow.jl") #support for forcing data from Delft3D-FM

include("era5.jl") #support for ERA5 reanalysis data

include("particles_core.jl") #support for background images downloaded with WMS

#unstructured_grid.jl
export Grid, Interpolator, add_grid!, interpolate, nodes_per_cell, winding_number, find_first_cell, get_values_by_cells!, find_first_cell, find_cells, create_node_tree!, dump, dump_subgrid, dump_array

#cartesian_grid.jl
export CartesianGrid, dump, in_bbox, find_index, find_index_and_weights, apply_index_and_weights, interpolate, CartesianXYTGrid, get_map_slice, update_cache, weights

#wms_client.jl
export WmsServer, get_map, plot_image

#dflow.jl
export load_nc_info, load_dflow_grid, load_nc_var, load_nc_map_slice, find_index, apply_index, get_times, get_reftime, as_DateTime, initialize_interpolation

#era5.jl
export EraData, load_map_slice, get_reftime, get_times, as_DateTime, initialize_interpolation

#particles_core.jl
export default_userdata, run_simulation, print_times, plot_maps_xy, plot_maps_xz, index

end #module
