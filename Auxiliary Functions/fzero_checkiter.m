function stop = fzero_checkiter(x,optimValues,state,maxevals)
    if optimValues.funccount >= maxevals
        stop = true;
    else
        stop = false;
    end
end