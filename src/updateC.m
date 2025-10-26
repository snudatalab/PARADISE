% Update the coefficients C for each slice (fixed basis Z).
%
% INPUT:
% Z: fixed basis matrices Z
% US: seasonal factors
% lambda_s, lambda_l: regularization constants for seasonality and L2
%
% OUTPUT
% C_result: updated factor matrices C
function C_result = updateC(Z, US, lambda_s, lambda_l)

    R = size(US, 2);
    K = length(US);


    parfor k = 1:K
        ZtZ = Z{k}' * Z{k};
        M = lambda_s * ZtZ + lambda_l * eye(size(ZtZ));
    
        % compute inverse using cholesky decompositionv
        [L, flag] = chol(M, 'lower');
        
        y = L \ (lambda_s * Z{k}' * US{k});
        C_result{k} = (L' \ y);
    end

end
