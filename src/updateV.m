% Update the shared factor V (column-wise, observed entries only).
%
% INPUT:
% X: input slices
% missing_ind_mat: sparse masks of missing entries (0=observed, 1=missing)
% UT, US, V, W: factor matrices UT, US, V, and W at the current iteration
% lambda_l: regularization constant for L2
%
% OUTPUT
% V_result: updated factor matrix V
function V_result = updateV(X, missing_ind_mat, UT,US, V, W, lambda_l)
    R = size(W, 2);
    J = size(V, 1);
    K = length(X);

    
    parfor k=1:K
        U{k} = UT{k} + US{k};
    end



    XUS = zeros(J, R);
    parfor k=1:K
        XUS = XUS +  (X{k}' * U{k} * diag(W(k,:)));
    end

    
    parfor j = 1:J
        utmp = zeros(K, R, R);
        for k = 1:K
            nonzero_ind = find(missing_ind_mat{k}(:,j) == 0);
            Utmp = U{k}(nonzero_ind, :);               
            utmp(k,:,:) = diag(W(k,:)) * (Utmp' * Utmp) * diag(W(k,:));  
        end
        
        % compute inverse using cholesky decomposition
        M = lambda_l * eye(R) + squeeze(sum(utmp, 1));
    
        [L, flag] = chol(M, 'lower');
        
        y = (L \ XUS(j, :)');
        V_result(j, :) = (L' \ y)';
    end
end