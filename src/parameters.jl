abstract type LineType end

@with_kw struct StFunction <: LineType
    τ::Float64 = 55.
    tstart::Float64 = 0.001
    tend::Float64 = 3.
    tN::Int64 = 425

    #These parameters are functions of the previous. Cannot specify explicitly.
    t::Vector{Float64} = ones(1)
    S::Function = (t) -> (sqrt.(t/(pi*τ)).*(exp.(-τ ./t) .- 1).+erf.(sqrt.(τ ./t)))
    function StFunction(τ,tstart,tend,tN,t,S)
        t=[LinRange( Float64(tstart), Float64(tend), Int64(tN));]
        new(τ,tstart,tend,tN,t,S)
    end
end

abstract type KMC <: LineType end

@with_kw struct KMC_HT <: KMC
    #These parameters do not change within an optimization run
    L::Int64 = 100 #Rod length
    z0::Int64 = 10 #Initial position
    nwalkers::Int64 = 50000 #Number of walkers

    #These parameters are redefined on each iteraction of the optimization
    #But they do not change for a given walk
    khop::Float64 = 0.4 #Hopping rate between lattice sites
    p::Float64 = 0.6 #Probability of site being occupied by an acceptor
    kHT::Float64 = 113. #Transfer rate to acceptor

    #These parameters are functions of the previous. Cannot specify explicitly.
    Γ::SArray{Tuple{3},Float64,1,3} = SVector{3,Float64}(1.,1.,1.) #Cumulative rate vector
    Γinv::SArray{Tuple{3},Float64,1,3} = 1.0 ./ Γ #Cumulative waiting time vector

    function KMC_HT(L,z0,nwalkers,khop,p,kHT,Γ,Γinv)
        Γ = SVector{3,Float64}(khop,2.0*khop,2.0*khop+kHT)
        Γinv = 1.0 ./ Γ
        new(L,z0,nwalkers,khop,p,kHT,Γ,Γinv)
    end
end

@with_kw struct KMC_RC <: KMC
    #These parameters do not change within an optimization run
    L::Int64 = 100 #Rod length
    z0::Int64 = 10 #Initial position
    nwalkers::Int64 = 50000 #Number of walkers

    #These parameters are redefined on each iteraction of the optimization
    #But they do not change for a given walk
    khop::Float64 = 0.4 #Hopping rate between lattice sites
    Γ::SArray{Tuple{2},Float64,1,2} = SVector{2,Float64}(1.,1.) #Cumulative rate vector
    Γinv::SArray{Tuple{2},Float64,1,2} = Γ #Cumulative waiting time vector

    function KMC_RC(L,z0,nwalkers,khop,Γ,Γinv)
        Γ = SVector{2,Float64}(khop,2.0*khop)
        Γinv = 1.0 ./ Γ
        new(L,z0,nwalkers,khop,Γ,Γinv)
    end
end

@with_kw struct Data <: LineType
    xdata::Vector{Float64} = ones(1)
    ydata::Vector{Float64} = ones(1)
    filename::String = "CdSODPA_Bulb.csv"
   #whererun::String = "local" # local or summit  ~/.julia/dev/DukovicSA two vs one dot
    function Data(xdata,ydata,filename) #See: https://stackoverflow.com/questions/45328831/what-on-earth-is-an-inner-constructor
        a = string("../",filename)
        xdata,ydata = load_ordered(joinpath(@__DIR__, a))
        new(xdata,ydata,filename)
    end
end