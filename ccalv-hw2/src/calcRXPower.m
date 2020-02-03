%Cameron J. Calv
function rxpower_in_dbm = calcRXPower(mobilePos, servingCellCenter, frequency ,function_request)

dist_real = real(mobilePos)-real(servingCellCenter);
dist_imag = imag(mobilePos)-imag(servingCellCenter);
distance = sqrt(dist_real^2+dist_imag^2);


if strcmp(function_request,'friis')
%         distance = distance (m)
%         txPower = 1 (mW) 
%         frequency = frequency (Hz)
%         txGain = 1
%         rxGain = 1
%         systemLosses = 1 (mW)
    rxpower_in_dbm = mw2dbm(friisFreeSpace(distance, 1, frequency, 1, 1, 1));
elseif strcmp(function_request, 'path_loss_exponent_3')
%         distance = distance (m)
%         reference_power_dbm = 0 (dBm)
%         reference_distance = 1 (m)
%         path_loss_exponent = 3
%         shadow_stdev = 8 (dBm)
%         shadowing_on = 0 (False)
    rxpower_in_dbm = pathLossExponent(distance, 0, 1, 3, 8, 0);
elseif strcmp(function_request, 'path_loss_exponent_4')
%         distance = distance (m)
%         reference_power_dbm = 0 (dBm)
%         reference_distance = 1 (m)
%         path_loss_exponent = 4
%         shadow_stdev = 8 (dBm)
%         shadowing_on = 0 (False)
    rxpower_in_dbm = pathLossExponent(distance, 0, 1, 4, 8, 0);
elseif strcmp(function_request, 'path_loss_exponent_3_shadowing')
%         distance = distance (m)
%         reference_power_dbm = 0 (dBm)
%         reference_distance = 1 (m)
%         path_loss_exponent = 3
%         shadow_stdev = 8 (dBm)
%         shadowing_on = 1 (True)
    rxpower_in_dbm = pathLossExponent(distance, 0, 1, 3, 8, 1);
elseif strcmp(function_request, 'path_loss_exponent_4_shadowing')
%         distance = distance (m)
%         reference_power_dbm = 0 (dBm)
%         reference_distance = 1 (m)
%         path_loss_exponent = 4
%         shadow_stdev = 8 (dBm)
%         shadowing_on = 1 (True)
    rxpower_in_dbm = pathLossExponent(distance, 0, 1, 4, 8, 1);
end

end