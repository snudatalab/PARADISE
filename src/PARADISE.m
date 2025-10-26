% -------------------------------------------------------------------------
% INPUTS
% X: input slices
% Z: fixed sinusoid/Fourier basis
% R: rank
% missing_ind: missing/observed index metadata per slice
% maxiter, conv: maximum #iterations and early-stop tolerance
% lambda_u, lambda_t, lambda_s, lambda_l: regularization constants
%
% OUTPUTS
% U: temporal factor (trend + seasonal)
% S: slice weights
% V: shared column factor
% fit_each: per-iteration NRE on observed entries
% times: average update time per iteration (sec)
% -------------------------------------------------------------------------
function [U, S, V, fit_each, times] = PARADISE(X, Z, R, missing_ind, maxiter, conv, lambda_u, lambda_t, lambda_s, lambda_l)

%% initialization
    ConvCrit = conv;
    flag = 0;
    fit_each = zeros(maxiter, 1);
    K = length(X);

    J = size(X{1}, 2);
    U = cell(K,1);
    H = rand(R,R);
    V = rand(J,R);
    W = rand(K,R);
    S = cell(K, 1);
    Q = cell(K, 1);
    UT = cell(K, 1);
    US = cell(K, 1);
    for i = 1:K
       S{i} = diag(rand(R,1)); 
       Q{i} = randn(size(X{i}, 1), R);
       UT{i} = randn(size(X{i}, 1), R);
       US{i} = randn(size(X{i}, 1), R);
       C{i} = randn(size(Z{i}, 2), R);
    end

    missing_ind_mat = cell(K,1);
    parfor k=1:K
        Ik = size(X{k}, 1);
        J = size(X{k}, 2);
        missing_ind_mat{k} = sparse(missing_ind{k}(:,2), missing_ind{k}(:,3), 1, Ik, J);
    end

    normX = 0;
    for k=1:length(X)
        normX = normX + norm(X{k}(find(missing_ind_mat{k} == 0)), "fro")^2;
    end

    recLoss = 0;
    iter = 0;
    times = 0;
    %% training
    for it=1:maxiter
        oldrec = recLoss;

        %% update factor matrices
        tic;    
        UT = updateUT(X, missing_ind_mat, UT, US, Q, W, V, H, lambda_u, lambda_t);
        US = updateUS(X, missing_ind_mat, UT, Q, W, V, H, Z, C, lambda_u, lambda_s);
        C = updateC(Z, US, lambda_s, lambda_l);
        W = updateW(X, missing_ind_mat, UT, US, V, lambda_l);
        V = updateV(X, missing_ind_mat, UT, US, V, W, lambda_l);
        [Q, H] = updateQH(H, UT, US);
        times = times + toc;

        parfor k = 1:K
            U{k} = UT{k} + US{k};
        end

        recLoss = 0;

        for k=1:K
            Ik = size(X{k}, 1);
            diff = X{k} - U{k} * diag(W(k,:)) * V';
            diff = diff(find(missing_ind_mat{k} == 0));
            recLoss = recLoss + norm(diff, "fro")^2;
        end
        NRE = recLoss / normX;
        fit_each(it) = NRE;
        
		if ((it > 1) &&  (abs(recLoss-oldrec)<oldrec*ConvCrit))
            flag = 1;
        else
            flag = 0;
        end

        fprintf('Iter %2d: oldrec = %7.3e, newrec = %7.3e, NRE = %.4f\n', ...
            it, oldrec, recLoss, fit_each(it));  
        iter = iter + 1;
        if flag == 1
            break
        end
        
    end
    
    
    parfor k=1:K
        S{k} = diag(W(k,:));
    end
    times = times / iter;

   
end