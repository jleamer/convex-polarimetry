function [theta, d] = dop(filename, sheet_name)
%UNTITLED4 Summary of this function goes here
%   Function to find the minimium value of epsilon that 
%   Read in the data from the table
data = table2array(readtable(filename, 'Sheet',sheet_name));
num_rows = size(data,1);
theta = data(:,1);
   
d = zeros(num_rows, 1);
for i = 1:num_rows
    s0 = data(i, 2) + data(i, 3);
    s1 = data(i, 2) - data(i, 3);
    s2 = 2 * data(i, 4);
    s3 = 2 * data(i, 5);
    
    d(i, 1) = sqrt(s1^2 + s2^2 + s3^2)/s0;
    
end

end

