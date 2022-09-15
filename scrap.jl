
############ My shit
@profview DukovicSA.run_optim(optimsystem,testsystem,P0,LB,UB,RUNTIME=1);
@profview fitKMC(optimsystem,testsystem.xdata,testsystem.ydata,P0);

|using BenchmarkTools
462.080 ms
@btime DukovicSA.fitKMC(optimsystem,testsystem.xdata);
@btime DukovicSA.run_optim(optimsystem,testsystem,P0,LB,UB,RUNTIME=1);

@btime ne(1)
@btime gt(1)

using Random
rng = Random.MersenneTwister()
a = (-1,1)
@btime rand(rng,a)
gol() = rand(rng,a)
@btime gol()

@btime rand(rng,Float64)
@btime fitKMC(optimsystem,testsystem.xdata,testsystem.ydata,P0);

## Test Interpolation
##Data outside bounds##
xdata = Vector(LinRange(-0.1,1.25,40))
xmc = sort(rand(50))
ymc = rand(50)
plot(xmc,zero(xmc),label="KMC")
plot!(xdata,ones(size(xdata)),label="Data")

itp = interpolate((xmc,) , ymc, Gridded(Linear()))
etp = extrapolate(itp, NaN)
yint = etp(xdata)

plot(xmc,ymc,label="KMC")
plot!(xdata,yint,label="Data")

intacceptor = DukovicSA.int_data(xmc,ymc,xdata)

##Data outside bounds##
xdata = Vector(LinRange(-0.1,1.25,70))
xmc = sort(rand(50))
ymc = rand(50)
plot(xmc,zero(xmc),label="KMC")
plot!(xdata,ones(size(xdata)),label="Data")

itp = interpolate((xmc,) , ymc, Gridded(Linear()))
etp = extrapolate(itp, NaN)
yint = etp(xdata)

plot(xmc,ymc,label="KMC")
plot!(xdata,yint,label="Data")

##KMC outside bounds##
xdata = Vector(LinRange(0.1,0.8,40))
xmc = sort(rand(50))
ymc = rand(50)
plot(xmc,zero(xmc),label="KMC")
plot!(xdata,ones(size(xdata)),label="Data")

itp = interpolate((xmc,) , ymc, Gridded(Linear()))
etp = extrapolate(itp, NaN)
yint = etp(xdata)

plot(xmc,ymc,label="KMC")
plot!(xdata,yint,label="Data")

##Find bounds
s = 0
[s+=1 for (index,value) in enumerate(yint) if !isnan(yint)] 

f = Iterators.filter(!isnan, yint)
foreach(println, f)

s= rand(length(yint))
sum([(s[index] - x)^2 for (index,x) in enumerate(yint) if !isnan.(x)])

function fancyitr(s,yint)
    sum([(s[index] - x)^2 for (index,x) in enumerate(yint) if !isnan.(x)])
end
@btime fancyitr(s,yint)

idx = Array{Bool}(undef,size(s))
function notsofancy(yint,idx,s)
    a = 0
    idx .= isnan.(yint)
    idx .= .!(idx)
    for i in eachindex(idx)
        if idx[i]
            a += (s[i] - yint[i])^2 
        end
    end
    return a
end
@btime notsofancy(yint,idx,s)





d = sort(Distributions.rand(optimsystem.distribution(),optimsystem.nwalkers))
d ./= maximum(d)
d .*= maximum(testsystem.xdata)
dinv  = DukovicSA.makehist(d)
reverse!(dinv)
plot(d,dinv,xaxis=:log)
plot!(xfit,yfit*P0[1],l=2,label="β=1")

density(d,label="pdf",l=1,marker=2)


data = testsystem.ydata
idx = sample(1:length(data),length(data))
resdatay = data[sort(idx)]
resdatax = testsystem.xdata[sort(idx)]
plot!(resdatax,resdatay,xaxis=:log)
