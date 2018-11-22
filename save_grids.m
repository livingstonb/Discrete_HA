clear

sgrid = cell(1,18);
i = 0;
for nx = [50 100 150 200 500 2000]
    for curv = [0.2 0.3 0.4]
        i = i + 1;
        sgrid{i} = linspace(0,1,nx)';
        sgrid{i} = sgrid{i} .^ (1 ./ curv);
        sgrid{i} = 1000 .* sgrid{i};
        
        T = table(sgrid{i});
        writetable(T,['/Users/Brian/Documents/midway2_output/discrete_time/11_21_18/grid',num2str(i),'.xls'],'WriteRowNames',true);
    end
end
