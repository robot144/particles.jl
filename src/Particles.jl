
module Particles

using NetCDF

include("unstructured_grid.jl") #support for fast indexing of unstructured grids

include("wms_client.jl") #support for background images downloaded with WMS

include("dflow.jl") #support for background images downloaded with WMS

include("particles_core.jl") #support for background images downloaded with WMS

export Grid, Interpolator, add_grid!, interpolate, nodes_per_cell, winding_number, find_first_cell, get_values_by_cells!, find_first_cell, find_cells, create_node_tree!, dump, dump_subgrid, dump_array

export WmsServer, get_map, plot_image

export load_nc_info, load_dflow_grid, load_nc_var, load_nc_map_slice, find_index, apply_index, get_times, get_reftime, as_DateTime, initialize_interpolation

export default_userdata, run_simulation, print_times, plot_maps_xy, plot_maps_xz, index

end #module
