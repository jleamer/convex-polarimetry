clear all; close all; clc
%file_list = ["HWP_hand_high.xlsx", "HWP_polarimeter_high.xlsx", "HWP_hand_low.xlsx", "HWP_polarimeter_low.xlsx"];
file_list = ["CP_35mW_hand.xlsx", "CP_350mW_hand.xlsx", "CP_1.21mW_polarimeter.xlsx", "CP_16.6mW_polarimeter.xlsx"];
sheet_names = ["calculated", "rho_mat", "rho_min", "rho_max"];

for i = 1:4
    computeRho(file_list(i))
    [theta, d_calc] = dop(file_list(i), sheet_names(1));
    [theta, d_rho] = dop(file_list(i), sheet_names(2));
    [theta, d_rhomin] = dop(file_list(i), sheet_names(3));
    [theta, d_rhomax] = dop(file_list(i), sheet_names(4));
    figure(i)
    plot(theta, d_calc, '-o', 'LineWidth',1.25)
    hold on
    plot(theta, d_rho, '-*', 'LineWidth',1.25)
    plot(theta, d_rhomin, '-x', 'LineWidth',1)
    plot(theta, d_rhomax, '-+', 'LineWidth',1)
    %title("DOP for " + file_list(i), 'Interpreter','none')
    xlabel("\fontsize{16} \theta")
    ax = gca;
    ax.FontSize = 12;
    %ylim([0.75, 1.05])
    xticks([0:20:180])
    ylabel("\fontsize{16} DOP")
    legend({'measured', 'corrected', 'min purity', 'max purity'}, 'Fontsize',13, 'Position',[0 0 1 1])
    hold off
end