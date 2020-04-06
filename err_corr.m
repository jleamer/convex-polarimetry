file_list = ["HWP_hand_high.xlsx", "HWP_polarimeter_high.xlsx", "HWP_hand_low.xlsx", "HWP_polarimeter_low.xlsx"];
sheet_names = ["calculated", "rho_mat"];

for i = 1:4
    computeRho(file_list(i))
    [theta, d_calc] = dop(file_list(i), sheet_names(1));
    [theta, d_rho] = dop(file_list(i), sheet_names(2));
    figure(i)
    plot(theta, d_calc)
    hold on
    plot(theta, d_rho)
    title("DOP for " + file_list(i))
    legend("calculated", "corrected")
    hold off
end

%L1-norm is called compressed sensing