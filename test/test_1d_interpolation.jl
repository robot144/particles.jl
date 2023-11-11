#
# testing cartesian_grid.jl 
#
using NetCDF

#
# low-level tests
#

function test1() #increasing x values and focus on edge cases
    x_values = [0.0, 1.0, 2.0, 3.0]
    y_values = [0.0, 1.0, 2.0]
    x1 = -0.1 # out of range
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x1, extrapolate=true,order=0) == 0.0       
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x1, extrapolate=true,order=1) == 0.0 # constant extrapolation   
    @test isnan(interpolation_linear_grid_edge_value_center(x_values, y_values, x1, extrapolate=false,order=0))         
    @test isnan(interpolation_linear_grid_edge_value_center(x_values, y_values, x1, extrapolate=false,order=1))
    x2 = 0.0  # at the edge
    #@test interpolation_linear_grid_edge_value_center(x_values, y_values, x2, extrapolate=false,order=1) ≈ 0.0 # edge now undefined
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x2, extrapolate=true,order=1) ≈ 0.0
    x3 = 0.25 # left of the first x_center
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x3, extrapolate=false,order=1) ≈ 0.0
    x4 = 0.5  # at the first x_center
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x4, extrapolate=false,order=1) ≈ 0.0
    x5 = 0.75 # right of first x_center 
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x5, extrapolate=false,order=0) ≈ 0.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x5, extrapolate=false,order=1) ≈ 0.25
    x6 = 1.25
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x6, extrapolate=false,order=0) ≈ 1.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x6, extrapolate=false,order=1) ≈ 0.75
    x7 = 1.75
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x7, extrapolate=false,order=0) ≈ 1.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x7, extrapolate=false,order=1) ≈ 1.25
    x8 = 2.25
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x8, extrapolate=false,order=0) ≈ 2.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x8, extrapolate=false,order=1) ≈ 1.75
    x9 = 2.75
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x9, extrapolate=false,order=0) ≈ 2.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x9, extrapolate=false,order=1) ≈ 2.0
    x10 = 3.0 # at the edge
    #@test interpolation_linear_grid_edge_value_center(x_values, y_values, x10, extrapolate=false,order=0) ≈ 2.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x10, extrapolate=true,order=1) ≈ 2.0  
    x11 = 3.1 # out of range
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x11, extrapolate=true,order=0) ≈ 2.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x11, extrapolate=true,order=1) ≈ 2.0 # constant extrapolation
    @test isnan(interpolation_linear_grid_edge_value_center(x_values, y_values, x11, extrapolate=false,order=0))
    @test isnan(interpolation_linear_grid_edge_value_center(x_values, y_values, x11, extrapolate=false,order=1))  
end

function test2() #decreasing x values and focus on edge cases
    x_values = [3.0, 2.0, 1.0, 0.0]
    y_values = [2.0, 1.0, 0.0] # reversed order, but same interpolation
    x1 = -0.1 # out of range
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x1, extrapolate=true,order=0) == 0.0       
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x1, extrapolate=true,order=1) == 0.0 # constant extrapolation   
    @test isnan(interpolation_linear_grid_edge_value_center(x_values, y_values, x1, extrapolate=false,order=0))         
    @test isnan(interpolation_linear_grid_edge_value_center(x_values, y_values, x1, extrapolate=false,order=1))
    x2 = 0.0  # at the edge
    #@test interpolation_linear_grid_edge_value_center(x_values, y_values, x2, extrapolate=false,order=1) ≈ 0.0 # edge now undefined
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x2, extrapolate=true,order=1) ≈ 0.0
    x3 = 0.25 # left of the first x_center
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x3, extrapolate=false,order=1) ≈ 0.0
    x4 = 0.5  # at the first x_center
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x4, extrapolate=false,order=1) ≈ 0.0
    x5 = 0.75 # right of first x_center 
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x5, extrapolate=false,order=0) ≈ 0.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x5, extrapolate=false,order=1) ≈ 0.25
    x6 = 1.25
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x6, extrapolate=false,order=0) ≈ 1.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x6, extrapolate=false,order=1) ≈ 0.75
    x7 = 1.75
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x7, extrapolate=false,order=0) ≈ 1.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x7, extrapolate=false,order=1) ≈ 1.25
    x8 = 2.25
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x8, extrapolate=false,order=0) ≈ 2.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x8, extrapolate=false,order=1) ≈ 1.75
    x9 = 2.75
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x9, extrapolate=false,order=0) ≈ 2.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x9, extrapolate=false,order=1) ≈ 2.0
    x10 = 3.0 # at the edge
    #@test interpolation_linear_grid_edge_value_center(x_values, y_values, x10, extrapolate=false,order=0) ≈ 2.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x10, extrapolate=true,order=1) ≈ 2.0  
    x11 = 3.1 # out of range
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x11, extrapolate=true,order=0) ≈ 2.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x11, extrapolate=true,order=1) ≈ 2.0 # constant extrapolation
    @test isnan(interpolation_linear_grid_edge_value_center(x_values, y_values, x11, extrapolate=false,order=0))
    @test isnan(interpolation_linear_grid_edge_value_center(x_values, y_values, x11, extrapolate=false,order=1))
end

function test3() #increasing x values and not equaly spaced
    x_values = [0.0, 1.0, 4.0, 10.0]
    y_values = [0.0, 1.0, 2.0]
    x1 = 2.0 
    # constant in 2nd cell -> 1.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x1, extrapolate=true,order=0) ≈ 1.0       
    # linear between centers at 0.5 and 2.5 with values 0.0 and 1.0 -> 0.75
    # 0.5 1.0  1.5 2.0  2.5 x
    # 0.0 0.25 0.5 0.75 1.0 y
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x1, extrapolate=true,order=1) ≈ 0.75 #   
    x2 = 8.0 # right of the last x_center 
    # constant in last cell -> 2.0
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x2, extrapolate=true,order=0) ≈ 2.0       
    # also constant after last center
    @test interpolation_linear_grid_edge_value_center(x_values, y_values, x2, extrapolate=true,order=1) ≈ 2.0 #   
end 

function test4() #increasing x values and focus on edge cases
    x_values = [0.0, 1.0, 2.0, 3.0]
    y_values = [0.0, 1.0, 2.0, 3.0]
    x1 = -0.1 # out of range
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x1, extrapolate=true,order=1) ≈ 0.0 # constant extrapolation   
    @test isnan(interpolation_linear_grid_edge_value_edge(x_values, y_values, x1, extrapolate=false,order=1))
    x2 = 0.0  # at the edge
    #@test interpolation_linear_grid_edge_value_edge(x_values, y_values, x2, extrapolate=false,order=1) ≈ 0.0 # edge now undefined
    @test isapprox(interpolation_linear_grid_edge_value_edge(x_values, y_values, x2, extrapolate=true,order=1),0.0,atol=1e-15)
    x3 = 0.25 # left of the first x_edge
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x3, extrapolate=false,order=1) ≈ 0.25
    x4 = 0.5  # at the first x_edge
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x4, extrapolate=false,order=1) ≈ 0.5
    x6 = 1.25
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x6, extrapolate=false,order=1) ≈ 1.25
    x7 = 1.75
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x7, extrapolate=false,order=1) ≈ 1.75
    x8 = 2.25
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x8, extrapolate=false,order=1) ≈ 2.25
    x9 = 2.75
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x9, extrapolate=false,order=1) ≈ 2.75
    x10 = 3.0 # at the edge
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x10, extrapolate=true,order=1) ≈ 3.0  
    x11 = 3.1 # out of range
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x11, extrapolate=true,order=1) ≈ 3.0 # constant extrapolation
    @test isnan(interpolation_linear_grid_edge_value_edge(x_values, y_values, x11, extrapolate=false,order=1))  
end

function test5() #decreasing x values and focus on edge cases
    x_values = [3.0, 2.0, 1.0, 0.0]
    y_values = [3.0, 2.0, 1.0, 0.0] # reversed order, but same interpolation
    x1 = -0.1 # out of range
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x1, extrapolate=true,order=1) ≈ 0.0 # constant extrapolation   
    @test isnan(interpolation_linear_grid_edge_value_edge(x_values, y_values, x1, extrapolate=false,order=1))
    x2 = 0.0  # at the edge
    @test isapprox(interpolation_linear_grid_edge_value_edge(x_values, y_values, x2, extrapolate=true,order=1),0.0,atol=1e-15)
    x3 = 0.25
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x3, extrapolate=false,order=1) ≈ 0.25
    x4 = 0.5 
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x4, extrapolate=false,order=1) ≈ 0.5
    x5 = 0.75  
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x5, extrapolate=false,order=1) ≈ 0.75
    x6 = 1.25
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x6, extrapolate=false,order=1) ≈ 1.25
    x7 = 1.75
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x7, extrapolate=false,order=1) ≈ 1.75
    x8 = 2.25
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x8, extrapolate=false,order=1) ≈ 2.25
    x9 = 2.75
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x9, extrapolate=false,order=1) ≈ 2.75
    x10 = 3.0 # at the edge
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x10, extrapolate=true,order=1) ≈ 3.0  
    x11 = 3.1 # out of range
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x11, extrapolate=true,order=1) ≈ 3.0 # constant extrapolation
    @test isnan(interpolation_linear_grid_edge_value_edge(x_values, y_values, x11, extrapolate=false,order=1))
end

function test6() #increasing x values and not equaly spaced
    x_values = [0.0, 1.0, 4.0, 10.0]
    y_values = [0.0, 1.0, 4.0, 10.0]
    x1 = 2.0   
    # linear between 1.0 and 3.0 -> 2.0
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x1, extrapolate=true,order=1) ≈ 2.0 #   
    x2 = 8.0 
    @test interpolation_linear_grid_edge_value_edge(x_values, y_values, x2, extrapolate=true,order=1) ≈ 8.0 #   
end 

test1()
test2()
test3()

test4()
test5()
test6()

