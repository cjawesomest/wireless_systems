% ulaSteering.m
%     Cameron J. Calv

function [steering_vector, spatial_spectrum] = ulaSteering(userPositions, cellCenter, frequency, M, d)
    % d given in fractions of wavelength
    wavelength = (3e8)/frequency;
    k = (2*pi)/wavelength;
    antenna_locations = linspace(cellCenter, cellCenter-j*(M)*(d*wavelength), M+1);    
    antenna_locations = antenna_locations(1:M);
    
    syms steering_vector(theta)
    steering_vector(theta) = [];
    for idx=1:M
       steering_vector(theta) = [steering_vector(theta) exp((-1)*j*k*(idx-1)*d*wavelength*cos(theta))];
    end
    
    sig_x_s = [];
    R_sum = zeros(M, M);
    for idx=1:numel(userPositions)
        sig_x = [];
        for a_idx=1:M
            sig_x = [sig_x ...
                dbm2mw(calcRXPower(antenna_locations(a_idx), userPositions(idx), frequency, 'path_loss_exponent_2.9'))...
                *exp((-1)*j*k*(a_idx-1)*d*wavelength*cos(atan((real(userPositions(idx))...
            -real(antenna_locations(a_idx)))/(imag(antenna_locations(a_idx))...
            -imag(userPositions(idx))))))]; 
        end
        R_sum = R_sum + transpose(sig_x)*(transpose(sig_x)');
        sig_x_s = [sig_x_s; sig_x];
    end
    R = R_sum / numel(userPositions);
    normR = R - min(R(:));
    normR = normR ./ max(normR(:));
    R = normR;
    
    
    steering_vector(theta) = transpose(steering_vector(theta));
    spatial_spectrum(theta) = ((steering_vector(theta)')*R*steering_vector(theta))/((steering_vector(theta)')*steering_vector(theta));
    spatial_spectrum(theta) =  (1/2)*((steering_vector(theta)')*ones(M, M)*steering_vector(theta))/((steering_vector(theta)')*steering_vector(theta)) - spatial_spectrum(theta);
%     spatial_spectrum(theta) = ((steering_vector(theta)')*R*steering_vector(theta))/((steering_vector(theta)')*steering_vector(theta));    
    
%     figure(2);
%     fplot(real(spatial_spectrum), [0 pi]);
%     spatial_spectrum(theta) = sqrt(real(spatial_spectrum(theta))^2+imag(spatial_spectrum(theta))^2);
    spatial_spectrum(theta) = real(spatial_spectrum(theta));

end
