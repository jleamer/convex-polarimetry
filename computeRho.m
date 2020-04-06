function computeRho(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%   Read in data from filename
    sheet_name = 'rho_mat';
    data = table2array(readtable(filename, 'Sheet','calculated'));
    rows = size(data,1);
    results = table();
    varNames = ["theta", "Jxx", "Jyy", "beta", "gamma", "trace_sq"];
    
%   Iterate over the data and solve for rho
    for i = 1:rows
        J = [[data(i,2), complex(data(i,4), data(i,5))]; ...
             [complex(data(i,4), -data(i,5)), data(i,3)]];
         
        cvx_begin
            variable rho(2,2) semidefinite complex
            minimize(norm(rho-J, 2))
            subject to
                trace(rho) == 1
        cvx_end
        
        temp = table(data(i,1), rho(1,1), rho(2,2), real(rho(1,2)), imag(rho(1,2)), trace(rho)^2, 'VariableNames',varNames);
        results = [results;temp];
    end
    
    writetable(results, filename, 'Sheet',sheet_name)
end

