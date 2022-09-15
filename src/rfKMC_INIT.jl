function initializeholes(sys)
    a = round(Int,sys.holedist[1] - 1.7320508075688772*sys.holedist[2])
    b = round(Int,sys.holedist[1] + 1.7320508075688772*sys.holedist[2])
    return a:b
end

function makehist!(eventtimes,yplot_PTZ)
    if length(eventtimes)!=0
        #Histogram is 1:(# of events); normalized by # of events
        yplot_PTZ .= [1.0/length(eventtimes):1.0/length(eventtimes):1.0;] #Create normalized # of transferred holes to PTZ
    else
        yplot_PTZ .= copy(eventtimes)
    end
    return nothing
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
