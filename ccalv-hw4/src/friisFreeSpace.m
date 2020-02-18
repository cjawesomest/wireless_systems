%Cameron J. Calv
function rxPower_in_mw = friisFreeSpace(distance, txPower, frequency, txGain, rxGain, systemLosses)
    rxPower_in_mw = txPower*((txGain+rxGain)/systemLosses)*(((3*10^6)/frequency)/(4*pi*distance))^2;
end