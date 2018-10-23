function params = parameters_experiment(runopts)

    freq = 1;
    
    params = MPCParams(freq,'test');
    params.nxlong = 500;
    params.Display = 1;
    
    if runopts.fast == 1
        params.set_fast();
    end
    
    params.index = 1;
end