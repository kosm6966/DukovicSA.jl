using DukovicSA
using Plots

## Load TestSystem, typically data
testsystem = Data(
    filename="220801_194K_after1ns.csv",
    logt_start = 0,
    logt_end = 6,
);

## Set up optimization, choose to fit S(t), KMC_RT, or KMC_HT to TestSystem

##Bulb
optimsystem = KMC_RC(
    L=101, #Length of NNR in units of sites
    holedist = (51., 14.), #Distribution parameters
    nwalkers = 50000, #Number of walkers
    khop = 140., #Hopping rate between lattice sites
);

# P0 = [A,khop,μ,σ]
P0 = [0.12,768.45,82.32,10.83];
LB = [0., 0., 0.5, 0.2886751345948129];
UB = [2., 1000.,optimsystem.L-1., 0.2886751345948129*(optimsystem.L-1.)];

yfit,xfit = fitKMC(optimsystem,testsystem.xdata,testsystem.ydata,P0);
plot(testsystem.xdata,testsystem.ydata,l=2,xaxis=:log,yaxis=:log,label = "data",legend=:topright,xlabel="ns",ylabel="arb. units")
plot!(xfit,yfit*P0[1],l=2,label="SA")

# Run SA
run_optim(optimsystem,testsystem,P0,LB,UB;RUNTIME=10);

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


testsystem = Data(
    filename="211001_CdSODPA_PTZ_h.csv",
    logt_start = 0,
    logt_end = 4,
);

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
    holedist = (1,10), #Distribution parameters
    nwalkers = 50000, #Number of walkers
    khop = 0.4, #Hopping rate between lattice sites
    p = 0.6, #Probability of site being occupied by an acceptor
    kHT = 113. #Transfer rate to acceptor
);


for m = 1:testsystem.logt_end
    a = testsystem.logt_start:testsystem.logt_end
    sumlog = 0
    for n = 1:length(testsystem.xdata)
        if testsystem.xdata[n] > 10. ^a[m] && testsystem.xdata[n] < 10. ^a[m+1]
            sumlog += 1
        end
    end
    # print(a[m],',',a[m+1],'\n')    
    print(string("time interval: 1E",a[m],"-1E",a[m+1],", n points:"), sumlog,'\n')
end

### Plot it
function fitKMC(sys,xdata,ydata,P0;label1="data",label2="KMC",xlabel="ns",ylabel="arb. units")
    sys = KMC_RC(sys,khop=P0[2],holedist=(P0[3],P0[4])) #Store inputs in structure
    d = DukovicSA.initializeholes(sys)
	if d[1] < 1 || d[2] > optimsystem.L
		return (zero(ydata).+1), xdata
	end
    tacceptor = zeros(Float64,sys.nwalkers)
	yacceptor = zeros(Float64,sys.nwalkers)
    yfit = similar(ydata)
	DukovicSA.runwalkers!(sys,tacceptor)
    DukovicSA.makehist!(tacceptor,yacceptor)
    DukovicSA.reversebulbonly!(sys,yacceptor)
    yacceptor ./= maximum(yacceptor)
    yacceptor .+= 1.
    DukovicSA.int_data!(tacceptor,yacceptor,xdata,yfit)
    xfit = copy(xdata)
    DukovicSA.cutxy!(yfit,xfit)
	plt = plot(xdata,ydata,l=2,xaxis=:log,yaxis=:log,label=label1,legend=:topright,xlabel=xlabel,ylabel=ylabel);
	plt = plot!(xfit,yfit*P0[1],l=2,label=label2);
	return yfit, xfit
end

