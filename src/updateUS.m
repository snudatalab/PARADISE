% Update the seasonal factor U^s for each slice (row-wise, observed entries only).
%
% INPUT:
% X: input slices
% missing_ind_mat: sparse masks of missing entries (0=observed, 1=missing)
% UT: current trend factors
% Q, H, W, V, Z, C: factor matrices Q, H, W, V, Z, and C at the current iteration
% lambda_u, lambda_s: regularization constants for uniqueness and seasonality
%
% OUTPUT
% U_result: updated factor matrices Uks
function U_result = updateUS(X, missing_ind_mat, UT, Q, W, V, H, Z, C, lambda_u, lambda_s)
    R = size(W, 2);
    J = size(V, 1);
    K = length(X);
    
   
    parfor k = 1:K
        Ik = size(X{k}, 1);
        Uk = (X{k} * V * diag(W(k ,:))) + lambda_u * (Q{k} * H - UT{k}) + lambda_s * Z{k} * C{k};
        for i = 1:Ik
            nonzero_ind = find(missing_ind_mat{k}(i,:) == 0);
            inter = diag(W(k,:)) * (V(nonzero_ind,:)' * V(nonzero_ind,:)) * diag(W(k,:));
            Uk(i, :) = Uk(i, :) - UT{k}(i, :) * inter;

            % compute inverse using cholesky decomposition
            M = (lambda_u+ lambda_s) * eye(R) + inter;
            [L, flag] = chol(M, 'lower');
            y = L \ Uk(i, :)';
            Uk(i, :) = (L' \ y)';
        end

        U_result{k} = Uk;
    end

end