"""

`int_data(xmodel,ymodel,xdata)` interpolates model data and returns two vectors.

For nonuniformly spaced xdata, must use 'Gridded' interpolation.These use a very similar syntax to BSplines, with the major exception being that one does not get to choose the grid representation (they are all OnGrid). There are only three modes for Gridded:

Gridded(Linear()) - a linear interpolation is applied between nodes
Gridded(Constant()) - nearest neighbor interpolation is used on the applied axis
NoInterp() ?

See: http://juliamath.github.io/Interpolations.jl/latest/control/

"""
function int_data(x,y,newx)
        # for i in eachindex(x)
        #         if x[i] < newx[1] || x[i] > newx[end]
        #                 x[i] = NaN
        #         end
        # end
        # cutxy!(x,y)
        # For nonuniform spacing, pass the x-values as a tuple since this generalizes to multi-dimensional coordinates
        etp = extrapolate(interpolate((x,) , y, Gridded(Linear())), NaN)
        yint = etp.(newx)
        return vec(yint)
end

function cutxy!(y1,y2)
        a = findall(isnan,y1)
        deleteat!(y1,a)
        deleteat!(y2,a)
end
