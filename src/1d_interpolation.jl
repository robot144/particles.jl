#
# 1D Interpolation
#
# Simplified interpolation routine for 1D data. This is used for the vertical interpolation of the 3D data. 
# There the z-coordinates can change over space and time.
# Other spwcial cases are the interpolation of the waterlevel and the bedlevel, 
# i.e. the minimum and maximum of the z-range.

function fix_between(value,lower_bound,upper_bound)
    # fix value to be between lower_bound and upper_bound
    # if value is NaN, return NaN
    # if value is lower than lower_bound, return lower_bound
    # if value is higher than upper_bound, return upper_bound
    # otherwise return value
    if isnan(value)
        return NaN
    elseif value < lower_bound
        return lower_bound
    elseif value > upper_bound
        return upper_bound
    else
        return value
    end
end

function interpolation_linear_grid_edge_value_center(x_values::Vector, y_values::Vector, x::Float64; extrapolate=false, order=0)
    # linear interpolation of y(x) for a grid with x-values at the edges of the cells
    # and y-values at the center of the cells
    # x_values: x-values at the edges of the cells (not necessarily equidistant, but increasing or decreasing)
    # y_values: y-values at the center of the cells
    # extrapolate : true or false
    # order: 0 or 1 for constant or linear interpolation
    # x: x-value to interpolate
    # returnreturn NaN when out of range in case extrapolate=false
    # 
    # example:
    # x_values = [0.0, 1.0, 2.0, 3.0]
    # y_values = [0.0, 1.0, 2.0]
    # x = 0.5
    # y = interpolation_linear_grid_edge_value_center(x_values, y_values, x)
    # y == 0.0
    # x=0.1 -> y=0.0
    # x=1.0 -> y=0.5
    # x=3.5 -> y=NaN or 2.0

    decreasing = x_values[end] < x_values[1]
    i = searchsortedfirst(x_values, x, rev=decreasing)

    if i == 1 # first consider out of range cases
        if extrapolate
            return y_values[1]
        else
            return NaN
        end
    elseif i == length(x_values)+1
        if extrapolate
            return y_values[end]
        else
            return NaN
        end
    else # within range
        x1 = x_values[i-1]
        x2 = x_values[i]
        this_value = y_values[i-1]
        if order == 0
            return this_value
        elseif order == 1
            this_xc = (x1+x2)/2
            if decreasing
                if x<=this_xc
                    next_i = min(i+1, length(x_values))
                    next_xc = (x_values[next_i]+x_values[next_i-1])/2
                else
                    next_i = max(i-1, 2)
                    next_xc = (x_values[next_i]+x_values[next_i-1])/2
                end
            else
                if x<=this_xc
                    next_i = max(i-1, 2)
                    next_xc = (x_values[next_i]+x_values[next_i-1])/2
                else
                    next_i = min(i+1, length(x_values))
                    next_xc = (x_values[next_i]+x_values[next_i-1])/2
                end
            end
            next_value = y_values[next_i-1]
            next_weight = fix_between(abs(this_xc-x)/(abs(next_xc-this_xc)+eps(1.0)),0.0,1.0)
            this_weight = 1.0-next_weight
            #println("this $(i) $(this_value) $(this_xc) $(this_weight) next $(next_i) $(next_value) $(next_xc) $(next_weight)")
            return this_value*this_weight + next_value*next_weight
        else
            error("order must be 0 or 1")
        end
    end
end

function interpolation_linear_grid_edge_value_edge(x_values::Vector, y_values::Vector, x::Float64; extrapolate=false, order=1)
    # linear interpolation of y(x) for a grid with x-values at the edges of the cells
    # and y-values at the edge of the cells too
    # x_values: x-values at the edges of the cells (not necessarily equidistant, but increasing or decreasing)
    # y_values: y-values at the edge of the cells
    # extrapolate : true or false
    # order: 1 for linear interpolation
    # x: x-value to interpolate
    # returnreturn NaN when out of range in case extrapolate=false
    # 
    # example:
    # x_values = [0.0, 1.0, 2.0, 3.0]
    # y_values = [0.0, 1.0, 2.0, 3.0]
    # x = 0.5
    # y = interpolation_linear_grid_edge_value_edge(x_values, y_values, x)
    # y == 0.5
    # x=0.1 -> y=0.1
    # x=1.0 -> y=1.0
    # x=3.5 -> y=NaN or 3.5

    decreasing = x_values[end] < x_values[1]
    i = searchsortedfirst(x_values, x, rev=decreasing)

    if i == 1 # first consider out of range cases
        if extrapolate
            return y_values[1]
        else
            return NaN
        end
    elseif i == length(x_values)+1
        if extrapolate
            return y_values[end]
        else
            return NaN
        end
    else # within range
        x1 = x_values[i-1]
        x2 = x_values[i]
        y1 = y_values[i-1]
        y2 = y_values[i]
        if order == 1
            weight2 = fix_between(abs(x-x1)/(abs(x2-x1)+eps(1.0)),0.0,1.0)
            weight1 = 1.0-weight2
            #println("this $(i) $(this_value) $(this_xc) $(this_weight) next $(next_i) $(next_value) $(next_xc) $(next_weight)")
            return y1*weight1 + y2*weight2
        else
            error("order must be 0 or 1")
        end
    end
end
