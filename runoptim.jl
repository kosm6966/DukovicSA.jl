using DukovicSA
using Plots
using Distributions 
# See Distributions.jl docs for more info
# https://juliastats.org/Distributions.jl/latest/univariate/
# Choose Distribution: Uniform(a, b),  Normal(μ, σ), Exponential(θ), Weibull(α, θ)
# Choose what params to fit: :none, :first, :second, :both


## Load TestSystem, typically data
testsystem = Data(
    filename="CdSODPA_Bulb.csv"
);

testsystem = Data(
    filename="211001_CdSODPA_PTZ_h.csv"
);

## Set up optimization, choose to fit S(t), KMC_RT, or KMC_HT to TestSystem

##Analytical
# P0 = [A,τ]
P0 = [1.1, 19.];
LB = [0.1, 0.];

optimsystem = StFunction(
    τ=55., # Tau
    # Make a time vector on a log scale
    logt_start=-2, # Start time at 10^logt_start
    logt_end=5.3, # End time at 10^logt_end
    tpoints=500 # Number of time points
);

# Plot S(t)
plot(testsystem.xdata,testsystem.ydata,l=2,xaxis=:log,label = "testsystem",legend=:topright,xlabel="ns",ylabel="arb. units")
plot!(optimsystem.t,optimsystem.S(optimsystem.t)*P0[1],l=2,label = "optimsystem")

##Acceptors
# P0 = [A,khop,p,KHT;params]
P0 = [1.1, 19., 0.6,   32., 1];
LB = [0.1,  0.,  0.,    0., 0];
UB = [2., 1000., 1., 1000., 10];

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

##Bulb
# P0 = [A,khop;params]
P0 = [1., 1.]
LB = [0.1, 0.]
UB = [2., 1000.]

optimsystem = KMC_RC(
    L=100, #Length of NNR in units of sites
    distribution = Weibull, #Distribution for initial holes
    params = (1,2), #Distribution parameters
    fitparams  = :none, #Choose to fit Distribution parameters
    nwalkers = 50000, #Number of walkers
    khop = 4., #Hopping rate between lattice sites
);

# Plot distribution
plt,d = randomdistribution(optimsystem);
display(plt)

# Plot KMC
optimsystem = KMC_RC(optimsystem,params = (1,2))
plt,yfit,xfit = fitKMC(optimsystem,testsystem.xdata,testsystem.ydata,P0[1]);
yfit4 = yfit ./ maximum(yfit)
xfit4 = copy(xfit)
# display(plt)

plot(testsystem.xdata,testsystem.ydata,l=2,xaxis=:log,label = "data",legend=:topright,xlabel="ns",ylabel="arb. units")
plot!(xfit1,yfit1*1.1,l=2,label="β=1")

# Run SA
run_optim(optimsystem,testsystem,P0,LB,UB,RUNTIME=10);

## Save fit
data = [xfit yfit]
filename="KMC_uniform"
date="22_08_19"
saveascsv(data,string("Data/",date,"/",filename,".txt"))

## Save plot
filename="KMC_uniform"
date="22_08_19"
format="svg"
savefig(plt,string("Data/",date,"/",filename,".",format))

# Random Params
makerandp0file(LB,UB,6)