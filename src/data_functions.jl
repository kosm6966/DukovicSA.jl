"""

`int_data(xmodel,ymodel,xdata)` interpolates model data and returns two vectors.

For nonuniformly spaced xdata, must use 'Gridded' interpolation.These use a very similar syntax to BSplines, with the major exception being that one does not get to choose the grid representation (they are all OnGrid). There are only three modes for Gridded:

Gridded(Linear()) - a linear interpolation is applied between nodes
Gridded(Constant()) - nearest neighbor interpolation is used on the applied axis
NoInterp() ?

See: http://juliamath.github.io/Interpolations.jl/latest/control/

"""
function int_data(xmc,ymc,xdata)
        # for i in eachindex(x)
        #         if x[i] < newx[1] || x[i] > newx[end]
        #                 x[i] = NaN
        #         end
        # end
        # cutxy!(x,y)
        # For nonuniform spacing, pass the x-values as a tuple since this generalizes to multi-dimensional coordinates
        # etp = extrapolate(interpolate((x,) , y, Gridded(Linear())), NaN)
        # etp = interpolate((x,) , y, Gridded(Linear()))
        # yint = etp.(newx)
        etp = extrapolate(interpolate((xmc,) , ymc, Gridded(Linear())), NaN)
        return vec(etp(xdata))
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
        saveascsv(out,"../p0rand.txt")
end

function randomdistribution(sys;xlabel="z",ylabel="P(z)")
        d = Distributions.rand(sys.distribution(),sys.nwalkers)
        if minimum(d) < 0.5
            d .+= 0.5 - minimum(d)
        end
        plt = density(d,label="pdf",l=2)
        plt[1][1][:y] ./= maximum(plt[1][1][:y])

        d = Int64.(ceil.(d))
        plt = histogram!(d,normalize=:probability,label="binned",xlabel=xlabel,ylabel=ylabel,axticks = 0:1:(maximum(d)+(d[2]-d[1])))
        a = plt[1][2][:y]
        b = maximum(filter(!isnan,a))
        # plt[1][2][:y] ./= b
        plt[1][1][:y] .*= b

        yaxis!([0,b + 0.1])

        a = length(string(eltype(sys.distribution)))
        b = string(sys.distribution)[1:a]
        title!(string(b," ",string(sys.distribution)[(length(b)+a+3):end]))
        display(plt)
        return plt, d
end

# """
#
# `load_xy(x,y,filename)` saves the csv file at "filename" which contains columns of data, labeled x, y.
#
# """
# function load_xy(filename)
#         datapath = joinpath(@__DIR__, filename)
#         df = CSV.read(filename,DataFrames.DataFrame)
#         return vec(df.x), vec(df.y)
# end

"""

`load_ordered(x,y,filename)` opens the csv file at "filename" which contains columns of data, labeled x, y.

"""
function load_ordered(filename)
        datapath = joinpath(@__DIR__, filename)
        A = DelimitedFiles.readdlm(datapath, ',',Float64) # delimiters: comma ',' tab '\t'
        return vec(A[:,1]), vec(A[:,2])
end

# """
#
# `savecsv(x,y,filename)` saves the csv file at "filename" which contains columns of data, labeled x, y.
#
# """
# function saveascsv(x,y,filename)
# datapath = joinpath(@__DIR__, filename)
#     open(datapath,"w") do f
#         write(f, string("x, y\n"))
#         writedlm(f,[x y],",")
#     end
# 	return nothing
# end