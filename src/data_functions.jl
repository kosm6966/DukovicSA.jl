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
        # etp = extrapolate(interpolate((x,) , y, Gridded(Linear())), NaN)
        etp = interpolate((x,) , y, Gridded(Linear()))
        yint = etp.(newx)
        return vec(yint)
end

function cutxy!(y1,y2)
        a = findall(isnan,y1)
        deleteat!(y1,a)
        deleteat!(y2,a)
end

function cutoutofbounds!(xdata,ydata,xfit)
        idx = falses(length(xdata))
        for i in eachindex(xdata)
                idx[i] = (xdata[i] < xfit[1])
        end
        if idx[1]; deleteat!(xdata,idx); deleteat!(ydata,idx) end
        idx = falses(length(xdata))
        for i in eachindex(xdata)
                idx[i] = (xdata[i] > xfit[end])
        end
        if idx[end]; deleteat!(xdata,idx); deleteat!(ydata,idx) end
end

function saveascsv(A,filename)
	datapath = joinpath(@__DIR__, filename)
    open(datapath,"w") do f
        writedlm(f,A,",")
	end
end

function makerandp0file(LB,UB,n)
	out = zeros(n,length(LB))
	for i in eachindex(LB)
		out[:,i] = LB[i] .+ (UB[i] - LB[i]) .*rand(n)
	end
        out .= round.(out,digits=2)
        saveascsv(out,"p0rand.txt")
end