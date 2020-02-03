function rxPower = friisFreeSpace(distance, txPower, frequency, txGain, rxGain, systemLosses)
    rxPower = txPower*((txGain+rxGain)/systemLosses)*(((3*10^6)/frequency)/(4*pi*distance))^2;
end