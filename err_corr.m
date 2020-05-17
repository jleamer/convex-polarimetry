%file_list = ["HWP_hand_high.xlsx", "HWP_polarimeter_high.xlsx", "HWP_hand_low.xlsx", "HWP_polarimeter_low.xlsx"];
file_list = ["CP_35mW_hand.xlsx", "CP_350mW_hand.xlsx", "CP_1.21mW_polarimeter.xlsx", "CP_16.6mW_polarimeter.xlsx"];
sheet_names = ["calculated", "rho_mat", "rho_min"];

for i = 1:4
    computeRho(file_list(i))
    [theta, d_calc] = dop(file_list(i), sheet_names(1));
    [theta, d_rho] = dop(file_list(i), sheet_names(2));
    %[theta, d_rhomin] = dop(file_list(i), sheet_names(3));
    figure(i)
    plot(theta, d_calc)
    hold on
    plot(theta, d_rho)
    %plot(theta, d_rhomin)
    title("DOP for " + file_list(i), 'Interpreter','none')
    xlabel("\theta")
    ylabel("DOP")
    legend("measured", "unknown")
    hold off
end

%L1-norm is called compressed sensing