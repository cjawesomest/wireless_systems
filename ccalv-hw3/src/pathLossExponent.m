%Cameron J. Calv
function rxpower_in_dbm = pathLossExponent(distance, reference_power_dbm, ...
    reference_distance, path_loss_exponent, shadow_stdev, shadowing_on)
    
    if(shadowing_on)
        rxpower_in_dbm = reference_power_dbm-...
            10*path_loss_exponent*log(distance/reference_distance)/log(10)+...
            normrnd(0, shadow_stdev);
    else
        rxpower_in_dbm = reference_power_dbm-...
            10*path_loss_exponent*log(distance/reference_distance)/log(10);
    end
end

