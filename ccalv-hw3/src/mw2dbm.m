%Cameron J. Calv
function power_in_dbm = mw2dbm(power_in_milliwatts)
    power_in_dbm = 10 * log(power_in_milliwatts/1)/log(10);
end
