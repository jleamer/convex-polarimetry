function computeRho(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%   Read in data from filename
    sheet_name = 'rho_mat';
    sheet_name2 = 'rho_min';
    sheet_name3 = 'rho_max';
    data = table2array(readtable(filename, 'Sheet','calculated'));
    rows = size(data,1);
    varTypes = ["double", "double", "double", "double", "double", "double"];
    varNames = ["theta", "Jxx", "Jyy", "beta", "gamma", "trace_sq"];
    results = table('Size',[rows,6], 'VariableTypes',varTypes, 'VariableNames',varNames);
    min_results = table('Size',[rows,6], 'VariableTypes',varTypes, 'VariableNames',varNames);
    max_results = table('Size',[rows,6], 'VariableTypes',varTypes, 'VariableNames',varNames);
    eps = .1;
    identity = [[1 0 ];[0 1]];
    
%   Iterate over the data and solve for rho
    for i = 1:rows
        J = [[data(i,2), complex(data(i,4), data(i,5))]; ...
             [complex(data(i,4), -data(i,5)), data(i,3)]];
         
        cvx_begin quiet
            variable rho(2,2) semidefinite complex
            minimize(norm(rho-J, 'fro'))
            subject to
                trace(rho) == 1;
        cvx_end
        
        cvx_begin quiet
            variable rhomin(2,2) semidefinite complex
            variable tau
            minimize tau
            subject to 
                norm(rhomin - tau*J - 0.5*(1-tau)*identity, 'fro') <= 0
                norm(rhomin - J, 'fro') <= eps
                trace(rhomin) == 1
        cvx_end
       
        cvx_begin quiet
            variable rhomax(2,2) semidefinite complex
            variable tau
            maximize tau
            subject to
                norm(rhomax - tau*J - 0.5*(1-tau)*identity, 'fro') <= 0
                norm(rhomax - J, 'fro') <= eps
                trace(rhomax) == 1
        cvx_end
        
        temp = table(data(i,1), rho(1,1), rho(2,2), real(rho(1,2)), imag(rho(1,2)), trace(rho)^2, 'VariableNames',varNames);
        temp2 = table(data(i,1), rhomin(1,1), rhomin(2,2), real(rhomin(1,2)), imag(rhomin(1,2)), trace(rhomin)^2, 'VariableNames',varNames);
        temp3 = table(data(i,1), rhomax(1,1), rhomax(2,2), real(rhomax(1,2)), imag(rhomax(1,2)), trace(rhomax)^2, 'VariableNames',varNames);
        results(i,:) = temp;
        min_results(i, :) = temp2;
        max_results(i, :) = temp3;
    end
    
    writetable(results, filename, 'Sheet',sheet_name)
    writetable(min_results, filename, 'Sheet',sheet_name2)
    writetable(max_results, filename, 'Sheet',sheet_name3)
end

