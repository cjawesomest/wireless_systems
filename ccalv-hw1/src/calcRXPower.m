function rxPower = calcRXPower(mobilePos, servingCellCenter, frequency ,functionHandle)

dist_real = real(mobilePos)-real(servingCellCenter);
dist_imag = imag(mobilePos)-imag(servingCellCenter);
distance = sqrt(dist_real^2+dist_imag^2);

if(isa(functionHandle, 'function_handle'))
    f_info = functions(functionHandle);
    if(f_info.function == 'friisFreeSpace')
        rxPower = friisFreeSpace(distance, 1, frequency, 1/2, 1/2, 1);
    end
end
end
