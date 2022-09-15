using Distributions
using StatsBase

function fitKMC(sys::KMC,xdata,ydata,A;label1="data",label2="KMC",xlabel="ns",ylabel="arb. units")
	tacceptor = runwalkers(sys)
	yacceptor = makehist(tacceptor)
	reversebulbonly(sys,yacceptor)
	yfit = int_data(tacceptor,yacceptor,xdata)
	xfit = copy(xdata)
	cutxy!(yfit,xfit)
	# ydata =  ydata .+ abs.(minimum(testsystem.ydata)) .+ 0.001
	plt = plot(xdata,ydata,l=2,xaxis=:log,label=label1,legend=:topright,xlabel=xlabel,ylabel=ylabel);
	plt = plot!(xfit,yfit*A,l=2,label=label2);
	return plt, yfit, xfit
end

function initializeholes(rng,sys)
    dist = sys.distribution(sys.params...)
    d = rand(rng,dist,sys.nwalkers)
    d ./= maximum(d)
    d .*= sys.L
    return Int64.(ceil.(d))
end

abstract type KMC <: LineType end

@with_kw struct KMC_HT <: KMC
    #These parameters do not change within an optimization run
    L::Int64 = 100 #Rod length
    distribution = Uniform #Distribution for initial holes
    params::Tuple{Real,Real} = (1,10) #Distribution parameters
    fitparams::Symbol = :none #Choose to fit Distribution parameters
    nwalkers::Int64 = 50000 #Number of walkers

    #These parameters are redefined on each iteraction of the optimization
    #But they do not change for a given walk
    khop::Float64 = 0.4 #Hopping rate between lattice sites
    p::Float64 = 0.6 #Probability of site being occupied by an acceptor
    kHT::Float64 = 113. #Transfer rate to acceptor

    #These parameters are functions of the previous. Cannot specify explicitly.
    Γ::SArray{Tuple{3},Float64,1,3} = SVector{3,Float64}(1.,1.,1.) #Cumulative rate vector
    Γinv::SArray{Tuple{3},Float64,1,3} = 1.0 ./ Γ #Cumulative waiting time vector

    function KMC_HT(L,distribution,params,fitparams,nwalkers,khop,p,kHT,Γ,Γinv)
        Γ = SVector{3,Float64}(khop,2.0*khop,2.0*khop+kHT)
        Γinv = 1.0 ./ Γ
        new(L,distribution,params,fitparams,nwalkers,khop,p,kHT,Γ,Γinv)
    end
end

@with_kw struct KMC_RC <: KMC
    #These parameters do not change within an optimization run
    L::Int64 = 100 #Rod length
    distribution = Uniform #Distribution for initial holes
    params::Tuple{Real, Real} = (1,10) #Distribution parameters
    fitparams::Symbol = :none #Choose to fit Distribution parameters
    nwalkers::Int64 = 50000 #Number of walkers

    #These parameters are redefined on each iteraction of the optimization
    #But they do not change for a given walk
    khop::Float64 = 0.4 #Hopping rate between lattice sites
    Γ::SArray{Tuple{2},Float64,1,2} = SVector{2,Float64}(1.,1.) #Cumulative rate vector
    Γinv::SArray{Tuple{2},Float64,1,2} = Γ #Cumulative waiting time vector

    function KMC_RC(L,distribution,params,fitparams,nwalkers,khop,Γ,Γinv)
        Γ = SVector{2,Float64}(khop,2.0*khop)
        Γinv = 1.0 ./ Γ
        new(L,distribution,params,fitparams,nwalkers,khop,Γ,Γinv)
    end
end

function optimsysinit(P0,sys::KMC_HT)
	sys = KMC_HT(sys,khop=P0[2],p=P0[3],kHT=P0[4]) #Store inputs in structure
	if sys.fitparams == :first
		sys = KMC_HT(sys,params=(P0[5],sys.params[2]))
	elseif sys.fitparams == :second
		sys = KMC_HT(sys,params=(sys.params[1],P0[5]))
	elseif sys.fitparams == :both
		sys = KMC_HT(sys,params=(P0[5],P0[6]))
	end
	return sys
end

function optimsysinit(P0,sys::KMC_RC)
	sys = KMC_RC(sys,khop=P0[2]) #Store inputs in structure
	if sys.fitparams == :first
		sys = KMC_RC(sys,params=(P0[3],sys.params[2]))
	elseif sys.fitparams == :second
		sys = KMC_RC(sys,params=(sys.params[1],P0[3]))
	elseif sys.fitparams == :both
		sys = KMC_RC(sys,params=(P0[3],P0[4]))
	end
	return sys
end


optimsystem = KMC_RC(
    L=100, #Length of NNR in units of sites
    distribution = Weibull, #Distribution for initial holes
    params = (1,2), #Distribution parameters
    fitparams  = :none, #Choose to fit Distribution parameters
    nwalkers = 50000, #Number of walkers
    khop = 4., #Hopping rate between lattice sites
);

optimsystem = KMC_HT(
    L=100, #Length of NNR in units of sites
    distribution = Uniform, #Distribution for initial holes
    params = (1,10), #Distribution parameters
    fitparams  = :second, #Choose to fit Distribution parameters
    nwalkers = 50000, #Number of walkers
    khop = 0.4, #Hopping rate between lattice sites
    p = 0.6, #Probability of site being occupied by an acceptor
    kHT = 113. #Transfer rate to acceptor
);



function randomdistribution(sys;xlabel="z",ylabel="P(z)")
    d = Distributions.rand(sys.distribution(),sys.nwalkers)
    d ./= maximum(d)
    d .*= sys.L
    plt = density(d,label="pdf",l=2)
    plt[1][1][:y] ./= maximum(plt[1][1][:y])

    d = Int64.(ceil.(d))
    plt = histogram!(d,normalize=:probability,xaxis=:log,yaxis=:log,label="binned",xlabel=xlabel,ylabel=ylabel,axticks = 0:1:(maximum(d)+(d[2]-d[1])))
    a = plt[1][2][:y]
    b = maximum(filter(!isnan,a))
    # plt[1][2][:y] ./= b
    plt[1][1][:y] .*= b

    yaxis!([0,b + 0.01])

    a = length(string(eltype(sys.distribution)))
    b = string(sys.distribution)[1:a]
    title!(string(b," ",string(sys.distribution)[(length(b)+a+3):end]))
    display(plt)
    return plt, d
end

# using Distributions 
# See Distributions.jl docs for more info
# https://juliastats.org/Distributions.jl/latest/univariate/
# Choose Distribution: Uniform(a, b),  Normal(μ, σ), Exponential(θ), Weibull(α, θ)
# Choose what params to fit: :none, :first, :second, :both
using StatsPlots
using Plots

# Plot distribution
plt,d = randomdistribution(optimsystem);
display(plt)

# Plot KMC
optimsystem = KMC_RC(optimsystem,params = (1,2))
plt,yfit,xfit = fitKMC(optimsystem,testsystem.xdata,testsystem.ydata,P0[1]);
display(plt)