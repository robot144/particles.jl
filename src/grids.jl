#
# Generic grids and data_sources
#

# data_sources are the things that contain fields that can be interpolated
abstract type BaseDataSource end

"""
function initialize_interpolation(source::BaseDataSource,varname::String,reftime::DateTime
    ,dummy::Float64,cache_direction::Symbol=:forwards)
Return a function of x,y,x,t that interpolates the field varname.
This is a generic function that should be implemented for each type of DataSource    
"""
function initialize_interpolation(source::BaseDataSource,varname::String,reftime::DateTime
    ,dummy::Float64,cache_direction::Symbol=:forwards)
    error("initialize_interpolation not implemented for $(typeof(source))")
end

"""
function get_grid(source::BaseDataSource)
    Return the grid that underlies the data_source
"""
function get_grid(source::BaseDataSource)
    error("get_grid not implemented for $(typeof(source))")
end

# Base type for grids, which contain a discrete representation of a field
# They are the underlying data-structures for interpolation
abstract type BaseGrid end
abstract type SpaceGrid <: BaseGrid end     #not time dependent
abstract type SpaceTimeGrid <: BaseGrid end
