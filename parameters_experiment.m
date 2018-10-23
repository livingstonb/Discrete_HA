function params = parameters_experiment(runopts)

    freq = 4;
    
    params(1) = MPCParams(freq,'test');
    params(1).nxlong = 500;
    params(1).Display = 1;
    
    if runopts.fast == 1
        params(1).set_fast();
    end
    
    params(1).index = 1;
end