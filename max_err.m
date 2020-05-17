cp_file_list = ["CP_35mW_hand.xlsx", "CP_400mW_hand.xlsx", "CP_1.21mW_polarimeter.xlsx", "CP_16.6mW_polarimeter.xlsx"];
lin_file_list = ["HWP_hand_high.xlsx", "HWP_polarimeter_high.xlsx", "HWP_hand_low.xlsx", "HWP_polarimeter_low.xlsx"];
sheet_names = ["calculated", "rho_mat", "rho_min"];
cp_err = [0, 0, 0, 0];
lp_err = [0 0 0 0];
for i=1:4
    [theta, dcp_calc] = dop(cp_file_list(i), sheet_names(1));
    [theta, dcp_rho] = dop(cp_file_list(i), sheet_names(2));
    
    [theta, dlp_calc]= dop(lin_file_list(i), sheet_names(1));
    [theta, dlp_rho]= dop(lin_file_list(i), sheet_names(2));
    
    cp_err(i) = max(abs(dcp_calc - dcp_rho));
    lp_err(i) = max(abs(dlp_calc - dlp_rho));
end