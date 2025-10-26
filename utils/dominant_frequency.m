% Estimate dominant frequency per slice via FFT (time axis), averaging magnitudes across features.
%
% INPUT:
% X: input slices
%
% OUTPUT
% results: cell of structs {index, period, length} per slice
function results = dominant_frequency(X)
    K = length(X);
    results = cell(K, 1);

    for k = 1:K
        Xk = X{k};                 % T_k x F
        [T_k, ~] = size(Xk);

        % FFT along time axis for each feature
        Xk_fft = fft(Xk);          
        Xk_mag = abs(Xk_fft);     
        avg_mag = mean(Xk_mag, 2); 

        % dominant frequency (excluding DC component)
        [~, idx] = max(avg_mag(2:floor(T_k/2)));
        dom_idx = idx + 1;

        period = round(T_k / dom_idx);

        results{k} = struct( ...
            'index', dom_idx, ...
            'period', period, ...
            'length', T_k ...
        );
    end
end