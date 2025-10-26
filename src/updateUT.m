% Update the trend factor U^t for each slice (row-wise, observed entries only).
%
% INPUT:
%   X: input slices
%   missing_ind_mat: sparse masks of missing entries (0=observed, 1=missing)
%   UT_old: previous trend factors
%   US, Q, H, W, V: factor matrices US, Q, H, W, and V at the current iteration
%   lambda_u, lambda_t: regularization constant for uniqueness and trend
%
% OUTPUT
%   U_result: updated factor matrices Ukt
function U_result = updateUT(X, missing_ind_mat, UT_old, US, Q, W, V, H, lambda_u, lambda_t)
    R = size(W, 2);
    J = size(V, 1);
    K = length(X);
    
   
    parfor k = 1:K
        Ik = size(X{k}, 1);
        Vspecific = zeros(Ik, R, R);
        
        U_pad = [zeros(1, R);UT_old{k};zeros(1, R)];
        Uk = X{k} * V * diag(W(k, :)) + lambda_u * (Q{k} * H - US{k}) + lambda_t * (U_pad(1:end-2,:) + U_pad(3:end,:));

        for i=1:size(X{k}, 1)
            nonzero_ind = find(missing_ind_mat{k}(i,:) == 0);
            inter = diag(W(k,:)) * (V(nonzero_ind,:)' * V(nonzero_ind,:)) * diag(W(k,:));
            Uk(i, :) = Uk(i, :) - US{k}(i, :) * inter;
            if (i == 1) || (i == size(X{k}, 1))
                M = (lambda_u + lambda_t) * eye(R) + inter;
            else
                M = (lambda_u + 2 * lambda_t) * eye(R) + inter;
            end     

            % compute inverse using cholesky decomposition
            [L, flag] = chol(M, 'lower');
            y = L \ Uk(i, :)';
            Uk(i, :) = (L' \ y)';
        end
        U_result{k} = Uk;
    end

end