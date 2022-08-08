using DukovicSA
using Plots

optimsystem = StFunction() # P0 = [A,τ]
optimsystem = KMC_HT() # P0 = [A,khop,p,KHT;z0]
optimsystem = KMC_RC(z0=5,khop=1,nwalkers=50000) # P0 = [A,khop;z0]

testsystem = Data(filename="CdSODPA_Bulb.csv")

intacceptor, ydatacopy, xdata = DukovicSA.fitKMC(optimsystem,testsystem.xdata,testsystem.ydata);
plot(testsystem.xdata,testsystem.ydata,l=2,xaxis=:log,label = "NewCode",legend=:topright,xlabel="ns",ylabel="arb. units")
plot!(xdata,intacceptor*P0[1])

P0 = [1.0, 19.]
LB = [0.1,0.]
UB = [2.,1000.]
run_optim(optimsystem,testsystem,P0,LB,UB,RUNTIME=10)

## Save function
saveascsv([tacceptor yacceptor],"./Data/21_11_08/FullSim_211001_CdSODPA_PTZ_h.txt")

# Random Params
makerandp0file(LB,UB,6)
