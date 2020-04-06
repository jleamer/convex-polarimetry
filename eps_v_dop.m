filename = 'HWP_hand_high.xlsx';
sheet_name =  'calculated';
data = table2array(readtable(filename, 'Sheet',sheet_name));

powers = -2:1:10;
num = size(powers, 2);
eps = zeros(1, num);
for i = 1:num
    eps(i) = 10^(powers(i));
end

min_d = [];
for i = 1:5
    [theta, temp] = min_dop(filename, sheet_name, eps(i));
    min_d = [min_d; temp];
end

rel_err = (abs(min_d(1,:) - min_d(end, :)))./min_d(1,:);
plot(theta, rel_err)
xlabel('Theta (deg)')
ylabel('Rel Err')
title("Relative Error in DOP b/w eps = 1e" + powers(1) + " and 1e" + powers(end))