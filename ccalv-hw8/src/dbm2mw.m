%Cameron J. Calv
function power_in_milliwatts = dbm2mw(power_in_dbm)
    power_in_milliwatts = 1 * 10^(power_in_dbm/10);
end