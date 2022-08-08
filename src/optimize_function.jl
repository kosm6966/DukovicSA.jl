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

z0 = FitZ0()
"""
function run_optim(optimsystem::LineType,testsystem::LineType,P0,LB,UB;RUNTIME=86400,VERB=3,NT=50,NS=20,RT=0.85)

	xtest,ytest = DukovicSA.initializetest(testsystem)

	println(string("\n optimsystem: \n",optimsystem))
	println(string("\n testsystem: \n",testsystem))
	
	# Call optimize function from Optim package
	objective_function(P) = fit(P,optimsystem::LineType,ytest,xtest)
	calloptim(objective_function,LB,UB,P0,VERB,RUNTIME,NT,NS,RT)
	println("\n run complete \n")

	return nothing
end

function fit(P0,optimsystem::KMC,ytest,xtest)
	optimsystem = optimsysinit(P0,optimsystem)
	intacceptor, yint, _ = fitKMC(optimsystem,xtest,ytest)
	res = residual_function(P0[1],intacceptor,yint)
	return res
end

function fit(P0,optimsystem::StFunction,ytest,xtest)
	optimsystem = StFunction(optimsystem,τ=P0[2]) #Store inputs in structure
	STvec = optimsystem.S(xtest)
	res = residual_function(P0[1],ytest,STvec)
	return res
end

function optimsysinit(P0,sys::KMC_HT)
	sys = KMC_HT(sys,khop=P0[2],p=P0[3],kHT=P0[4]) #Store inputs in structure
	if length(P0) > 4 sys = KMC_HT(sys,z0=P0[5]) end#Store inputs in structure
	return sys
end

function optimsysinit(P0,sys::KMC_RC)
	sys = KMC_RC(sys,khop=P0[2]) #Store inputs in structure
	if length(P0) > 2 sys = KMC_HT(sys,z0=P0[3]) end#Store inputs in structure
	return sys
end

function residual_function(A,ytest,yfit)
	yfit .*= A
	yfit .-= ytest
	yfit .= yfit.^2
	return sum(yfit)/length(yfit)
end

function fitKMC(sys::KMC,xdata,ydata)
		tacceptor = DukovicSA.runwalkers(sys)
		yacceptor = DukovicSA.makehist(tacceptor)
		reversebulbonly(sys,yacceptor)
		ydatacopy = copy(ydata)
		xdatacopy = copy(xdata)
		DukovicSA.cutoutofbounds!(xdatacopy,ydatacopy,tacceptor)
		intacceptor = DukovicSA.int_data(tacceptor,yacceptor,xdatacopy)
		return intacceptor, ydatacopy, xdatacopy
end

function reversebulbonly(sys::KMC_RC,yvec)
	reverse!(yvec)

	return nothing
end

function reversebulbonly(sys::KMC_HT,yvec)

	return nothing
end

function initializetest(sys::KMC)
	tacceptor = runwalkersRC(sys)
	yacceptor = makehist(tacceptor)
	reversebulbonly(sys,yacceptor)

	return taccepter, yacceptor
end

function initializetest(sys::StFunction)
	return sys.t, sys.S(sys.t)
end

function initializetest(sys::Data)
	return copy(sys.xdata), copy(sys.ydata)
end

function calloptim(objective_function,LB,UB,P0,VERB,RUNTIME,NT,NS,RT)
	optimize(objective_function, LB, UB, P0,
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
end