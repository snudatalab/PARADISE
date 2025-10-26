delete(gcp('nocreate'))
numW = 20;
N = maxNumCompThreads(numW);
p = parpool('threads');

%% Load data
load_train_path = sprintf('./data/PEMS_train.mat');
load(sprintf(load_train_path, 'data'));
X = data;
load_test_path = sprintf('./data/PEMS_test.mat');
load(sprintf(load_test_path), 'data');
X_test = data;

K = length(X);

%% Prepare missing indices (metadata of observed entries) --------------
missing_ind = cell(1,K);
for k=1:K
    missing_ind{k} = X_test{k}(:, 1:3);  
end  

for k=1:length(X_test)
    cast(X_test{k}(:,1:3), 'int32');
    X_test{k}(:,1:3) = X_test{k}(:,1:3);
    nan_mask = isnan(X_test{k}(:, 4));
    X_test{k}(nan_mask, :) = [];
end

frequency = dominant_frequency(X);
P_list = cellfun(@(r) r.period, frequency);   
Z = fourier_basis(X, P_list, 5, 2);

%% Set hyperparameters
conv = 0.001; maxiter = 100;
R = 10; lambda_u = 0.01; lambda_l = 0.01; lambda_s = 0.01; lambda_t = 10;

%% train and evaluate model
[U, S, V, fit_each, times] = PARADISE(X, Z, R, missing_ind, maxiter, conv, lambda_u, lambda_t, lambda_s, lambda_l);
NRE_res = NRE(X_test, U, S, V, length(X));
fprintf('Proposed method fitness : %4f, time per iteration : %4f\n', NRE_res, times);  
