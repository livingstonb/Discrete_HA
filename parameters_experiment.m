function params = parameters_experiment(Fast)

    freq = 1;
    
    params = MPCParams(freq,'test');
    params.nxlong = 500;
    params.Display = 1;
    
    if Fast == 1
        params.set_fast();
    end

end