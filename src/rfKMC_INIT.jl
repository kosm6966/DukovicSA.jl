function initializeholes(rng,sys)
    dist = sys.distribution(sys.params...)
    d = rand(rng,dist,sys.nwalkers)
    d ./= maximum(d)*1/sys.L
    if minimum(d) < 0.5
        d .+= 0.5 - minimum(d)
    end
    return Int64.(ceil.(d))
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
