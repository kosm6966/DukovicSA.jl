@with_kw struct SystemParameters{T1<:Int64,T2<:Float64,T3<:SArray{Tuple{3},Float64,1,3}}
    #These parameters do not change within an optimization run
    L::T1 = 100 #Rod length
    z0::T1 = 3 #Initial position
    nwalkers::T1 = 20000 #Number of walkers

    #These parameters are redefined on each iteraction of the optimization
    #But they do not change for a given walk
    p::T2 = 0.6 #Probability of site being occupied by an acceptor
    khop::T2 = 12. #Hopping rate between lattice sites
    kHT::T2 = 113. #Transfer rate to acceptor
    Γ::T3 = SVector{3,Float64}(khop,2.0*khop,2.0*khop+kHT) #Cumulative rate vector
    Γinv::T3 = 1.0 ./ Γ #Cumulative waiting time vector
end

function makehist(eventtimes)
    if length(eventtimes)!=0
        #Histogram is 1:(# of events); normalized by # of events
        yplot_PTZ = collect(1.0/length(eventtimes):1.0/length(eventtimes):1.0) #Create normalized # of transferred holes to PTZ
    else
        yplot_PTZ = copy(eventtimes)
    end
    return yplot_PTZ
end

function gammawhat2(Γ1,uQ)
    choose = 1
    if uQ > Γ1; choose = 2; end
    return choose
end

function gammawhat3(Γ1,Γ2,uQ)
    choose = 1
    if uQ > Γ1; choose = 2; end
    if uQ > Γ2; choose = 3; end
    return choose
end
