"""

`run_optim(P0,LB,UB,Z0,NT,NS,RT,RUNTIME,VERB)` runs the Simulated Annealing
program to optimize the parameters given in P0 by recursively minimizing χ².
All inputs are required. Results are printed to screen.

P0: Array{Float64,1}. Initial Values for Parameters. [N,khop,KHT]
LB: Array{Float64,1}. Lower Bounds on Parameters. [lb1,lb2,lb3]
UB: Array{Float64,1}. Upper Bounds on Parameters. [ub1,ub2,ub3]
NT: Int64. Reduce temperature every NT*NS*dim(P0) evaluations.
NS: Int64. Adjust bounds every NS*dim(P0) evaluations.
RT: Float64. Cooling Rate. When temp changes, new temp is T_new = RT * T_old
RUNTIME: Int64. Terminate run after RUNTIME seconds.
VERB: Int64. If 0 - don't print to screen. If 3 - print everything to screen.
Z0: Int64. Either static z0 or z0max for uniform distribution

"""
function run_optim(P0,LB,UB,Z0,NT,NS,RT,RUNTIME,VERB)
#Uncomment if fitting KMC to data
#Data is stored in folder ~/.julia/dev/DukovicSA/
	lengthrod=100
	numwalkers=50000
	# xdata,ydata = load_ordered(joinpath(@__DIR__, "../211001_CdSODPA_PTZ_h.csv"))
	xdata,ydata = DukovicSA.load_ordered(joinpath(@__DIR__, "../CdSODPA_Bulb.csv"))

#Uncomment if fitting KMC_RC to S(t)
	# lengthrod=100
	# numwalkers=50000
	# tau = 55.
	# S = (t) -> (sqrt.(t/(pi*tau)).*(exp.(-tau ./t) .- 1).+erf.(sqrt.(tau ./t)))

#Uncomment if fitting S(t) to data
#Data is stored in folder ~/.julia/dev/DukovicSA/
	# S = (tau,t) -> (sqrt.(t/(pi*tau)).*(exp.(-tau ./t) .- 1).+erf.(sqrt.(tau ./t)))
	# xdata,ydata = DukovicSA.load_ordered(joinpath(@__DIR__, "../CdSODPA_Bulb.csv"))

#Uncomment this section if fitting S(t) to KMC_RC with static z0.
	# lengthrod=100
	# numwalkers=50000
	# kHOP=0.4
	# S = (tau,t) -> (sqrt.(t/(pi*tau)).*(exp.(-tau ./t) .- 1).+erf.(sqrt.(tau ./t)))
	# sysRC = DukovicSA.SystemParameters(L=lengthrod,z0=Z0,nwalkers=numwalkers,khop=kHOP) #Store parameters in structure
	# tacceptor = DukovicSA.runwalkersRC(sysRC) #For static z0
	# yacceptor = DukovicSA.makehist(tacceptor) #Make histogram of 1/t
	# reverse!(yacceptor) #Reverse for RC data

#Uncomment this section if fitting S(t) to KMC_RC with uniform z0.
	# lengthrod=100
	# numwalkers=50000
	# kHOP=0.4
	# S = (tau,t) -> (sqrt.(t/(pi*tau)).*(exp.(-tau ./t) .- 1).+erf.(sqrt.(tau ./t)))
	# sysRC = DukovicSA.SystemParameters(L=lengthrod,z0=Z0,nwalkers=numwalkers,khop=kHOP) #Store parameters in structure
	# tacceptor = DukovicSA.runwalkersRC_uniform(sysRC) #For uniform z0
	# yacceptor = DukovicSA.makehist(tacceptor) #Make histogram of 1/t
	# reverse!(yacceptor) #Reverse for RC data

	function objective_function(P0)
	#Uncomment if fitting KMC_RC to data with static z0. P0 = [khop,A].
		#sysRC = DukovicSA.SystemParameters(L=lengthrod,z0=Z0,nwalkers=numwalkers,khop=P0[1]) #Store inputs in structure
		#tacceptor = DukovicSA.runwalkersRC(sysRC)
		#yacceptor = DukovicSA.makehist(tacceptor)
		#reverse!(yacceptor)
		#intacceptor = DukovicSA.int_data(tacceptor,yacceptor,xdata)
		#yint = copy(ydata)
		#cutxy!(intacceptor,yint)
		#residual =  sum((P0[2]*intacceptor .- yint).^2)/length(yint)

	#Uncomment if fitting KMC_RC to data with uniform z0. P0 = [khop,A].
		sysRC = DukovicSA.SystemParameters(L=lengthrod,z0=Int64(round(P0[3])),nwalkers=numwalkers,khop=P0[1]) #Store inputs in structure
		tacceptor = DukovicSA.runwalkersRCv_uniform(sysRC)
		yacceptor = DukovicSA.makehist(tacceptor)
		reverse!(yacceptor)
		intacceptor = DukovicSA.int_data(tacceptor,yacceptor,xdata)
		yint = copy(ydata)
		cutxy!(intacceptor,yint)
		residual =  sum((P0[2]*intacceptor .- yint).^2)/length(yint)

	#Uncomment if fitting KMC_RC to S(t) with static z0. P0 = [khop,A].
		#sysRC = DukovicSA.SystemParameters(L=lengthrod,z0=Z0,nwalkers=numwalkers,khop=P0[1]) #Store inputs in structure
		#tacceptor = DukovicSA.runwalkersRC(sysRC)
		#yacceptor = DukovicSA.makehist(tacceptor)
		#reverse!(yacceptor)
		#residual =  sum((S(tacceptor) .- P0[2]*yacceptor).^2)/length(yacceptor)

	#Uncomment if fitting KMC_RC to S(t) with uniform z0. P0 = [khop,A].
		#sysRC = DukovicSA.SystemParameters(L=lengthrod,z0=Z0,nwalkers=numwalkers,khop=P0[1]) #Store inputs in structure
		#tacceptor = DukovicSA.runwalkersRC_uniform(sysRC)
		#yacceptor = DukovicSA.makehist(tacceptor)
		#reverse!(yacceptor)
		#residual =  sum((S(55.,tacceptor) .- P0[2]*yacceptor).^2)/length(yacceptor)

	#Uncomment if fitting S(t) to KMC_RC. P0 = [tau,A].
		#residual =  sum((P0[2]*S(P0[1],tacceptor) .- yacceptor).^2)/length(yacceptor)

	#Uncomment if fitting S(t) to data. P0 = [tau,A].
		#residual =  sum((P0[2]*S(P0[1],tacceptor) .- ydata).^2)/length(ydata)

	#Uncomment if fitting KMC_HT to data with static z0. P0 = [p,khop,kHT,A].
		#sysHT = SystemParameters(L=lengthrod,z0=Z0,nwalkers=numwalkers,p=P0[1],khop=P0[2],kHT=P0[3]) #Store inputs in structure
		#tacceptor = DukovicSA.runwalkersHT(sysHT)
		#yacceptor = DukovicSA.makehist(tacceptor)
		#intacceptor = DukovicSA.int_data(tacceptor,yacceptor,xdata)
		#yint = copy(ydata)
		#cutxy!(intacceptor,yint)
		#residual =  sum((P0[4]*intacceptor .- yint).^2)/length(yint)

	#Uncomment if fitting KMC_HT to data with uniform z0. P0 = [p,khop,kHT,A].
		#sysHT = SystemParameters(L=lengthrod,z0=Z0,nwalkers=numwalkers,p=P0[1],khop=P0[2],kHT=P0[3]) #Store inputs in structure
		#tacceptor = DukovicSA.runwalkersHT_uniform(sysHT)
		#yacceptor = DukovicSA.makehist(tacceptor)
		#intacceptor = DukovicSA.int_data(tacceptor,yacceptor,xdata)
		#yint = copy(ydata)
		#cutxy!(intacceptor,yint)
		#residual =  sum((P0[4]*intacceptor .- yint).^2)/length(yint)

        return residual
	end
# Call optimize function from Optim package
	res = optimize(objective_function, LB, UB, P0,
		SAMIN(;
			nt = NT,     # 10; reduce temperature every NT*NS*dim(P0) evaluations
			ns = NS,     # 20; adjust bounds every NS*dim(P0) evaluations
			rt = RT,     # 0.85; geometric temperature reduction factor: when temp changes, new temp is T=RT*T0
			neps = 5,  # number of previous best values the final result is compared to
			f_tol = 1e-4, # 1e-3, the required tolerance level for function value comparisons
			x_tol = 1e-4, # 1e-2, the required tolerance level for x
			coverage_ok = false, # if false, increase temperature until initial parameter space is covered
			verbosity = VERB), # scalar: 0, 1, 2 or 3 (default = 0).
		Optim.Options(
			iterations = 10^6,
	  		time_limit = RUNTIME, # seconds; 6 days: 518400
			show_trace = true,
			store_trace = false,
			extended_trace = true))

	println("\n run complete \n")

	return nothing
end
