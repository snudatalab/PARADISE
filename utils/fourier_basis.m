% Build fixed Fourier basis per slice using per-slice period P_k.
%
% INPUT:
% X: input slices
% P_list: per-slice periods (cell or numeric vector)
% K_freq, step: number of harmonics and spacing (default step=1)
%
% OUTPUT
% Z: Fourier basis matrices per slice (sin/cos for selected harmonics)
function Z = fourier_basis(X, P_list, K_freq, step)
    if nargin < 4
        step = 1; % default: consecutive harmonics
    end
    
    K = length(X);
    Z = cell(1, K);
    freq_list = 1:step:(step * K_freq);
    num_freqs = length(freq_list);

    for k = 1:K
        Xk = X{k};
        I_k = size(Xk, 1);
        t = (1:I_k)';

        if iscell(P_list)
            P_k = P_list{k};
        else
            P_k = P_list(k);
        end

        omega_k = 2 * pi / P_k;

        Z_k = zeros(I_k, 2 * num_freqs);
        for i = 1:num_freqs
            f = freq_list(i);
            Z_k(:, i) = sin(f * omega_k * t);
            Z_k(:, num_freqs + i) = cos(f * omega_k * t);
        end

        Z{k} = Z_k;
    end
end