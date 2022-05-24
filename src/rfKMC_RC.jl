#Notes from literature. Disregard.
# tau = 2-240
# A - Normalization
# f(t) = A*(sqrt(t/(Pi*tau))*(exp(-tau/t)-1)+erf(sqrt(tau/t)))
# tau - 21
# L - Use mean value of NNR lengths measured. Nature, 2016
#       Weighted by fraction with bulb: .32*22 + .59*37 + .65*68 /3 = 24
#       Non-weighted 22 + 37 + 68 /3 = 42
# SpecialFunctions.erf


function runwalkersRC(sys) #Returns the hitting times for nwalkers.
    #Pre-allocate arrays
    rng = Random.MersenneTwister() #Initalize seed
    tt = zeros(Float64,sys.nwalkers) #Store waiting times

    #Take n walks
    for n in 1:sys.nwalkers
        t = walkRC(sys,rng,sys.z0) #Static z0
        tt[n] = t
    end

    #Return n waiting times
    return sort!(tt)
end

function runwalkersRC_uniform(sys) #Returns the hitting times for nwalkers.
    #Pre-allocate arrays
    rng = Random.MersenneTwister() #Initalize seed
    tt = zeros(Float64,sys.nwalkers) #Store waiting times
    uniformz0 = rand(rng,1:sys.z0,sys.nwalkers) #'Look up table' for initial positions

    #Take n walks
    for n in 1:sys.nwalkers
        t = walkRC(sys,rng,uniformz0[n]) #Uniform z0
        tt[n] = t
    end

    #Return n waiting times
    return sort!(tt)
end

function walkRC(sys,rng,z0)
    t = 0. #Initialize walk time
    z = z0 #Initialize position
    RC = true #Initialize electron-hole recombination flag; Run til false
    go = event -> event==1 ? -1 : +1 #Go left(-) if event==1; if event==2
    bc = x -> x==(sys.L+1) ? (sys.L-1) : x #If event takes z from L to L+1 go to L-1

    while RC #Stop when it recombines
        choose = gammawhat2(sys.Γ[1],rand(rng,Float64)*sys.Γ[2]) #Choose event; Position of Γtotal=2 (- or +)
        z += go(choose) #Move walker. Γ1<uΓ<Γ2
        z = bc(z) #Impose boundary condition
        t += -log(rand(rng,Float64))*sys.Γinv[2] #Wait time
        RC *= (z!=0) #Did hole reencounter electron? If !=0 then do nothing (mult by true); If !(!=0) i.e. ==0 then RC=(true)*false=false; break loop
    end
    return t
end
