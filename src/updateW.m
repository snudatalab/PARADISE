% Update the slice weights W (per slice, observed entries only).
%
% INPUT:
% X: input slices
% missing_ind_mat: sparse masks of missing entries (0=observed, 1=missing)
% UT, US, V: factor matrices UT, US, and V at the current iteration
% lambda_l: regularization constant for L2
%
% OUTPUT
% W_result: updated slice weights W
function W_result = updateW(X, missing_ind_mat, UT, US, V, lambda_l)
    R = size(V, 2);
    J = size(V, 1);
    K = length(X);
    
    parfor k = 1:K
        U{k} = UT{k} + US{k};
    end

    parfor k = 1:K
        Ik = size(U{k}, 1);
        [nnz_ind_row, nnz_ind_col] = find(missing_ind_mat{k} == 0);
    
        A = V(nnz_ind_col,:) .* U{k}(nnz_ind_row,:);  
    
        M = lambda_l * eye(R) + A' * A;               
    
        % compute inverse using cholesky decomposition
        [L, flag] = chol(M, 'lower');
    
        vecVU = zeros(1, R);
        for r = 1:R
            vecVU(r) = U{k}(:,r)' * X{k} * V(:,r);
        end
    
        y = (L \ vecVU');
        W_result(k, :) = (L' \ y)';
    end
end