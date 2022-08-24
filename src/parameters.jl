abstract type LineType end

@with_kw struct StFunction <: LineType
    τ::Float64 = 55.
    logt_start::Real = -2
    logt_end::Real = 5.3
    tpoints::Int64 = 500

    #These parameters are functions of the previous. Cannot specify explicitly.
    t::Vector{Float64} = ones(1)
    S::Function = (t) -> (sqrt.(t/(pi*τ)).*(exp.(-τ ./t) .- 1).+erf.(sqrt.(τ ./t)))
    function StFunction(τ,logt_start,logt_end,tpoints,t,S)
        t = 10 .^(range(logt_start,stop=logt_end,length=tpoints))
        new(τ,logt_start,logt_end,tpoints,t,S)
    end
end

abstract type KMC <: LineType end

@with_kw struct KMC_HT <: KMC
    #These parameters do not change within an optimization run
    L::Int64 = 100 #Rod length
    distribution = Uniform #Distribution for initial holes
    params::Tuple{Real, Real} = (1,10) #Distribution parameters
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

@with_kw struct Data <: LineType
    xdata::Vector{Float64} = ones(1)
    ydata::Vector{Float64} = ones(1)
    filename::String = "CdSODPA_Bulb.csv"
    function Data(xdata,ydata,filename) #See: https://stackoverflow.com/questions/45328831/what-on-earth-is-an-inner-constructor
        a = string("../",filename)
        xdata,ydata = load_ordered(joinpath(@__DIR__, a))
        new(xdata,ydata,filename)
    end
end