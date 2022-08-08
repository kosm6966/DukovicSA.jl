using DukovicSA
using Plots

## Run rfKMC_RC
using SpecialFunctions
xdata,ydata = DukovicSA.load_ordered(joinpath(@__DIR__, "./CdSODPA_Bulb.csv"))

P0 = [0.4,1.]
Z0 = 10
sysRC = DukovicSA.SystemParameters(L=100,z0=Z0,nwalkers=500,khop=P0[1]) #Store inputs in structure
tacceptor = DukovicSA.runwalkersRC(sysRC)
# tacceptor = DukovicSA.runwalkersRC_uniform(sysRC)
yacceptor = DukovicSA.makehist(tacceptor); reverse!(yacceptor)

fgr
A = 0.86
plot(xdata,ydata,l=2,xaxis=:log,label = "CdSODPA_Bulb",legend=:topright,xlabel="s",ylabel="arb. units")
plot!(tacceptor,yacceptor,l=2,xaxis=:log,label = string("rfKMC_RC khop=",P0[1]),legend=:topright,xlabel="ns",ylabel="arb. units")
# savefig(string("Data/21_12_02/uniform_torfKMC_z0_",z0[iter],".svg"))


##Cut Interpolated KMC w/ example
function cutxy!(y1,y2,x)
        a = findall(isnan,y1)
        deleteat!(y1,a)
        deleteat!(y2,a)
		deleteat!(x,a)
end

yint = DukovicSA.int_data(tacceptor,yacceptor,xdata)
x = copy(xdata)
y = copy(ydata)
cutxy!(yint,y,x)
plot(x,y,l=2,xaxis=:log,label = "CdSODPA_Bulb",legend=:topright,xlabel="ns",ylabel="arb. units")
plot!(x,yint,l=2,xaxis=:log,label = "rfKMC_RC",legend=:topright,xlabel="s",ylabel="arb. units")

##Run SA
LB=[0.,0.];UB=[1000,2.];
Z0=5
#run_optim(P0,LB,UB,Z0,NT,NS,RT,RUNTIME,VERB)
DukovicSA.run_optim(Vector(P0),LB,UB,Z0,50,20,.85,10,3)

## Save function
using DelimitedFiles
function saveascsv(A,filename)
	datapath = joinpath(@__DIR__, filename)
    open(datapath,"w") do f
        writedlm(f,A,",")
	end
end
saveascsv([tacceptor yacceptor],"./Data/21_11_08/FullSim_211001_CdSODPA_PTZ_h.txt")

# Random Params
function randb(LB,UB,n)
	out = zeros(n,length(LB))
	for i in eachindex(LB)
		out[:,i] = LB[i] .+ (UB[i] - LB[i]) .*rand(n)
	end
	return out
end
#p0 = [z0,khop,kHT]
LB = [0.,0.]
UB = [1000.,2.]
n = 6
p0rand = round.(randb(LB,UB,n),digits=2)
saveascsv(p0rand,"p0rand.txt")
