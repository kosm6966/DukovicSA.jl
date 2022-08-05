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
