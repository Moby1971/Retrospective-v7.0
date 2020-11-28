function traject = array_trajectory(parameters)

% k-space trajectory

nr_pe = max(parameters.gp_var_mul)-min(parameters.gp_var_mul)+1;
nr_pes = parameters.NO_VIEWS;
pe_array = parameters.gp_var_mul;

traject=pe_array(1:nr_pes) + nr_pe/2 + 1;

end