function runwalkersHT(sys) #Returns the hitting times for nwalkers.
    rng = Random.MersenneTwister() #Initalize seed

    #Pre-allocate arrays
    rod = rand(rng,Float64,sys.L,sys.nwalkers) #Random numbers 'look up table' for placing acceptors
    idx = zeros(Int64,sys.L+2) #Used to store indices of acceptors
    isacceptor =  Vector{Bool}(undef,sys.L) #Stores 1 (true) if acceptors at position

    tt = zeros(Float64,sys.nwalkers) #Stores waiting times
    event = Vector{Bool}(undef,sys.nwalkers) #Stores type of event

    #Walk
    for n in 1:sys.nwalkers
        acceptors = init_acceptors!(idx,isacceptor,rod,sys,n) #Initialize acceptor positions
        t,whahahappened = walkHT(acceptors,sys,rng,sys.z0) #Do one walk
        tt[n] = t #Store waiting time of walk 'n'
        event[n] = whahahappened #If 1, HT; if 0 RC (recombine with e at bulb)
    end

    #Only take HT events and sort
    return sort!(tt[event])
end

function runwalkersHT_uniform(sys) #Returns the hitting times for nwalkers.
    rng = Random.MersenneTwister() #Initalize seed

    #Pre-allocate arrays
    rod = rand(rng,Float64,sys.L,sys.nwalkers) #Random numbers 'look up table' for placing acceptors
    uniformz0 = rand(rng,1:sys.z0,sys.nwalkers) #'Look up table' for initial positions
    idx = zeros(Int64,sys.L+2) #Used to store indices of acceptors
    isacceptor =  Vector{Bool}(undef,sys.L) #Stores 1 (true) if acceptors at position

    tt = zeros(Float64,sys.nwalkers) #Stores waiting times
    event = Vector{Bool}(undef,sys.nwalkers) #Stores type of event

    #Walk
    for n in 1:sys.nwalkers
        acceptors = init_acceptors!(idx,isacceptor,rod,sys,n) #Initialize acceptor positions
        t,whahahappened = walkHT(acceptors,sys,rng,uniformz0[n]) #Do one walk
        tt[n] = t #Store waiting time of walk 'n'
        event[n] = whahahappened #If 1, HT; if 0 RC (recombine with e at bulb)
    end

    #Only take HT events and sort
    return sort!(tt[event])
end

function init_acceptors!(A,isacceptor,rod,sys,iter)
    for i = 1:sys.L
        isacceptor[i] = (rod[i,iter] <= sys.p) #Acceptor at site 'i' if rand num < prob
    end
    #Acceptor positions in A. But need to include 0 and L+1. Position 0 is
    # the bulb. Position L+1 invokes the reflective boundary condition.
    # A = [0 (acceptor indices ∈ (1,L) ) L+1]
    # A[1] = 0
    # A[end] = L+1
    n = sum(isacceptor) + 2
    A[n] = sys.L+1
    A[2:n-1] = (1:sys.L)[isacceptor]
    return view(A,1:n) #Return view (pointer) to pre-allocated array on each iteration.
end

function walkHT(acceptors,sys,rng,z0)
    t = 0. #Initialize walk time
    z = z0 #Initialize position
    HT = true #Initialize hole transfer flag; Run until false
    RC = true #Initialize electron-hole recombination flag; Run until false
    go = event -> event==1 ? -1 : +1 #Go left(-) if event==1; if event==2 or event==3 move right(+) (event==3 kills loop anyhow)
    bc = x -> x==(sys.L+1) ? (sys.L-1) : x #If event takes z from L to L+1 go to L-1
    dt = numevents -> -log(rand(rng,Float64))*sys.Γinv[numevents] #Wait time; for non-acceptor site, numevents=2 for acceptor site, numevents=3

    #Walk left and right until you get to a site where an event can occur.
    bin = searchsorted(acceptors,z) #What is z between (acceptors (1-L), bulb(0), reflection point(L+1) )? Returns boundary indices.
    while HT && RC #Stop when walker recombines or transfers
        while z > acceptors[bin.stop] && z < acceptors[bin.start]
            choose = gammawhat2(sys.Γ[1],rand(rng,Float64)*sys.Γ[2]) #Choose event; Position of Γtotal=2 (- or +)
            z += go(choose) #Move walker. Γ1<uΓ<Γ3
            z = bc(z) #Impose boundary condition
            t += dt(2) #Add time; number of events = 2.
        end
        #If while loop breaks, at acceptor or bulb?
        RC *= (z!=0) #Did hole reencounter electron? Evaluate before next move incase the last move placed hole at bulb.
        choose = gammawhat3(sys.Γ[1],sys.Γ[2],rand(rng,Float64)*sys.Γ[3]) #Choose event; Position of Γtotal=3 (-, + or HT)
        z += go(choose) #Move walker
        z = bc(z) #Impose boundary condition
        t += dt(3)*(z!=0) #Add time if not at bulb; number of events = 3.
        HT *= (choose!=3) #Did HT occur? If !=3 do nothing (mult by true); If !(!=3) i.e. ==3 then HT=(true)*false=false; break loop
        RC *= (z!=0) #Did hole reencounter electron? If !=0 then do nothing (mult by true); If !(!=0) i.e. ==0 then RC=(true)*false=false; break loop
        if z < acceptors[bin.stop] || z > acceptors[bin.start] #If it moved into a new bin, what is the new bin?
            bin = bin .+ go(choose) #Update bin; only consider moving to adjacent bins
        end
        if length(bin)!=0
            bin = searchsorted(acceptors,z) #Finding new bin in above if statement is faster (I think) but doesn't work if z0 is on an acceptor site
        end
    end
    return t, !HT
end
