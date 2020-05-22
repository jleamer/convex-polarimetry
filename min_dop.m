function [theta, d] = min_dop(filename, sheetname, eps_guess)
%UNTITLED Summary of this function goes here
%   Function to find minimum value of epsilon as an error bound

%First read in the data
data = table2array(readtable(filename, 'Sheet',sheetname));
rows = size(data, 1);
varTypes = ["double", "double", "double", "double", "double"];
varNames = ["theta", "Jxx", "Jyy", "beta", "gamma"];
results = table('Size',size(data), 'VariableTypes',varTypes, 'VariableNames',varNames);
d = zeros(1, rows);

%Create the coherency matrix from data and then compute rho
for i = 1:rows
        J = [[data(i,2), complex(data(i,4), data(i,5))]; ...
             [complex(data(i,4), -data(i,5)), data(i,3)]];
         
        cvx_begin
            variables tau rho(2,2) semidefinite complex
            minimize(norm(rho-J, 2))
            subject to
                trace(rho) == 1
                norm(rho-J) <= eps_guess
        cvx_end
  
    
        s0 = rho(1, 1) + rho(2, 2);
        s1 = rho(1, 1) - rho(2, 2);
        s2 = 2 * real(rho(1, 2));
        s3 = 2 * imag(rho(1, 2));
    
        d(1, i) = sqrt(s1^2 + s2^2 + s3^2)/s0;
    
end

theta = data(:, 1);
end
