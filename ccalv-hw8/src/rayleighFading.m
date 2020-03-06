function [time, envelope, total_time] = rayleighFading(max_doppler_shift, num_samples)
    fm = max_doppler_shift;
    N = num_samples;
    delta_f = (2*fm)/(N-1);
    waveform_time = 1/(delta_f);

    noise_freqs = linspace(-fm, fm, N+1);

    % In_Phase
    noise_sources = [];
    for i=1:(N/2)
        noise_sources = [noise_sources (normrnd(0, 1)+j*normrnd(0, 1))];
    end
    noise_sources = [noise_sources 1 flip(conj(noise_sources))];

    doppler_spectrum = noise_sources.*sqrt(1.5./(pi.*fm.*sqrt(1-((noise_freqs./fm).^2))));
    last_slope = (doppler_spectrum(N)-doppler_spectrum(N-1))/(noise_freqs(N)-noise_freqs(N-1));
    doppler_spectrum(N+1) = doppler_spectrum(N)+last_slope;
    doppler_spectrum(1) = conj(doppler_spectrum(N+1));

    in_phase = ifft(ifftshift(doppler_spectrum));

    % Quadrature
    noise_sources = [];
    for i=1:(N/2)
        noise_sources = [noise_sources (normrnd(0, 1)+j*normrnd(0, 1))];
    end
    noise_sources = [noise_sources 1 flip(conj(noise_sources))];

    doppler_spectrum = noise_sources.*sqrt(1.5./(pi.*fm.*sqrt(1-((noise_freqs./fm).^2))));
    last_slope = (doppler_spectrum(N)-doppler_spectrum(N-1))/(noise_freqs(N)-noise_freqs(N-1));
    doppler_spectrum(N+1) = doppler_spectrum(N)+last_slope;
    doppler_spectrum(1) = conj(doppler_spectrum(N+1));

    quadrature = ifft(ifftshift(doppler_spectrum));
    
    time_series = sqrt((in_phase.^2) + (quadrature.^2));
    
    total_time = waveform_time;
    time = linspace(0, waveform_time, numel(time_series));
    envelope = 10*log10(time_series);
end