% Update Q and H (orthogonal alignment + averaged update).
%
% INPUT:
% H: factor matrix H
% UT, US: factor matrices UT and US at the current iteration
%
% OUTPUT
% Q_update, H_update
function [Q_update, H_update] = updateQH(H, UT, US)

    K = length(UT);
    R = size(H, 1);
    Q_update = cell(K, 1);
    H_tmp = zeros(K, R, R);
    H_update = zeros(R,R);

  

    parfor k=1:K
        U{k} = UT{k} + US{k};
        [Z, ~, P] = svd(U{k} * H', 'econ');
        Q_update{k} = Z * P';

        H_update = H_update +  Q_update{k}' * U{k};
    end

    H_update = H_update ./ K;

end

