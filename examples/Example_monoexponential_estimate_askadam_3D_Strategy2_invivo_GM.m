% Example_monoexponential_estimate_askadam_3D_Strategy2_invivo_GM
%
% Description:
%   Demonstrates fitting of a multi-flip-angle multi-echo monoexponential
%   decay model to 3D (simulated) in-vivo grey-matter-like data using the
%   askAdam solver (Strategy 2, memory-efficient masking approach).
%
%   Model:  S(alpha_k, TE) = S0_k .* exp(-TE .* (R2s_hat + dR2s_dalpha .* alpha_k))
%
%   S0 is flip-angle-dependent: each flip angle k has an independent S0_k
%   with no parametric constraint on how S0 varies with alpha.
%   R2s_hat and dR2s_dalpha are shared (common) across all flip angles,
%   i.e. the R2* model parameters are the same regardless of which flip
%   angle dataset is being fitted.
%
%   Priors for S0_k: uniform box priors [0, 2] are used for all flip angles.
%   This encodes only positivity and an upper bound consistent with a
%   unit-normalised acquisition without assuming any flip-angle dependence.
%
% Kwok-Shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 7 May 2025
% Date last modified:
%
addpath(genpath('../../gacelle/'))
clear

%% -----------------------------------------------------------------------
%  1. Simulate data
% -----------------------------------------------------------------------
% reproducibility
seed = 5438973; rng(seed); gpurng(seed);

% Flip angles [rad]  (e.g. 5, 10, 20, 30 degrees — typical mGRE vFA range)
alpha_deg = [5, 10, 20, 30];
alpha     = deg2rad(alpha_deg);
Na        = numel(alpha);

% Echo times [s]
Nt = 15;
t  = linspace(0, 40e-3, Nt);

% Volume size and SNR
Nx  = 21;
Ny  = 21;
Nz  = 21;
SNR = 100;

% Spherical brain mask
mask = strel('sphere', 10); mask = mask.Neighborhood;

% --- Ground-truth parameters ---
% R2* baseline [1/s] (GM: ~20-35 1/s)
R2s_hat_GT      = 30 + 5 * randn(Nx, Ny, Nz);
% Flip-angle sensitivity of R2* [1/s/rad]
dR2s_dalpha_GT  = 5  + 1 * randn(Nx, Ny, Nz);

% S0_k for each flip angle — no parametric form assumed.
% Here we choose values consistent with a simplified Ernst-angle model to
% generate realistic GT, but the fitting does NOT assume this form.
%   S0(alpha) ~ M0 * sin(alpha) * (1-E1) / (1 - cos(alpha)*E1)
%   with M0=1, T1~1.3 s (GM), TR~50 ms => E1=exp(-50e-3/1.3)
M0  = 1 + randn(Nx, Ny, Nz) * 0.2;    % proton-density-weighted baseline
E1  = exp(-50e-3 / 1.3);               % longitudinal recovery (fixed TR)

S0_GT = cell(1, Na);
for k = 1:Na
    S0_GT{k} = M0 .* sin(alpha(k)) .* (1-E1) ./ (1 - cos(alpha(k))*E1);
end

% --- Build GT parameter struct and generate forward signal ---
pars_GT.R2s_hat     = R2s_hat_GT;
pars_GT.dR2s_dalpha = dR2s_dalpha_GT;
for k = 1:Na
    pars_GT.(sprintf('S0_%d', k)) = S0_GT{k};
end

S_GT = Example_monoexponential_FWD_askadam_3D_Strategy2_GM(pars_GT, t, alpha, mask);

% Add Gaussian noise
noise_std = mean(cellfun(@(x) mean(abs(x(mask > 0))), S0_GT)) / SNR;
y         = S_GT + noise_std * randn(size(S_GT));

%% -----------------------------------------------------------------------
%  2. Set up askAdam fitting
% -----------------------------------------------------------------------

% --- Parameter names ---
% S0_k (one per flip angle), then shared R2* parameters
S0_param_names  = arrayfun(@(k) sprintf('S0_%d', k), 1:Na, 'UniformOutput', false);
modelParams     = [S0_param_names, {'R2s_hat', 'dR2s_dalpha'}];

% --- Priors (box constraints) ---
% S0_k:        [0, 2]    — positive, normalised signal; same for all FA
%                          (no FA dependence encoded in the prior)
% R2s_hat:     [0, 100]  — GM R2* baseline [1/s]
% dR2s_dalpha: [-10, 10] — flip-angle sensitivity [1/s/rad]
lb = [zeros(1, Na),  0,  -10];   % lower bounds (one per parameter)
ub = [2*ones(1, Na), 100,  10];  % upper bounds (one per parameter)

% --- Starting points ---
pars0 = struct();
for k = 1:Na
    pars0.(sprintf('S0_%d', k)) = 0.5 * ones(Nx, Ny, Nz) + 0.1 * randn(Nx, Ny, Nz);
end
pars0.R2s_hat     = 20 + 10 * randn(Nx, Ny, Nz);
pars0.dR2s_dalpha = zeros(Nx, Ny, Nz);

% --- Algorithm settings ---
fitting                     = [];
fitting.modelParams         = modelParams;
fitting.lb                  = lb;
fitting.ub                  = ub;
fitting.iteration           = 4000;
fitting.initialLearnRate    = 0.001;
fitting.lossFunction        = 'l1';
fitting.tol                 = 1e-4;
fitting.convergenceValue    = 1e-8;
fitting.convergenceWindow   = 20;
fitting.isDisplay           = false;
fitting.isOptimiseMemory    = true;

% Forward model
modelFWD = @Example_monoexponential_FWD_askadam_3D_Strategy2_GM;
weights  = [];

%% -----------------------------------------------------------------------
%  3. Run estimation
% -----------------------------------------------------------------------
askadam_obj = askadam;
out = askadam_obj.optimisation(y, mask, weights, pars0, fitting, modelFWD, t, alpha, mask);

%% -----------------------------------------------------------------------
%  4. Plot results
% -----------------------------------------------------------------------
sl = 11;   % display slice

% --- Scatter plots: GT vs estimated ---
figure('Name', 'S0 per flip angle — scatter');
tiledlayout(2, Na);
for k = 1:Na
    fname = sprintf('S0_%d', k);
    nexttile;
    scatter(S0_GT{k}(mask > 0), pars0.(fname)(mask > 0), 4, 'filled'); hold on;
    scatter(S0_GT{k}(mask > 0), out.final.(fname)(mask > 0), 4, 'filled');
    refline(1, 0);
    xlabel('GT'); ylabel(sprintf('S0_%d (alpha=%.0f°)', k, alpha_deg(k)));
    legend('Start', 'Fitted', 'Location', 'best');
end
nexttile;
scatter(R2s_hat_GT(mask > 0),     pars0.R2s_hat(mask > 0),     4, 'filled'); hold on;
scatter(R2s_hat_GT(mask > 0),     out.final.R2s_hat(mask > 0),     4, 'filled');
refline(1, 0); xlabel('GT'); ylabel('R2s\_hat [1/s]');
legend('Start', 'Fitted', 'Location', 'best');

nexttile;
scatter(dR2s_dalpha_GT(mask > 0), pars0.dR2s_dalpha(mask > 0), 4, 'filled'); hold on;
scatter(dR2s_dalpha_GT(mask > 0), out.final.dR2s_dalpha(mask > 0), 4, 'filled');
refline(1, 0); xlabel('GT'); ylabel('dR2s\_dalpha [1/s/rad]');
legend('Start', 'Fitted', 'Location', 'best');

% --- Image maps ---
figure('Name', 'Parameter maps');
Nrows = Na + 2;
tiledlayout(3, Nrows);
for k = 1:Na
    fname = sprintf('S0_%d', k);
    nexttile; imshow(S0_GT{k}(:,:,sl) .* mask(:,:,sl), [0 0.2]);
    title(sprintf('S0\\_%d GT (%.0f°)', k, alpha_deg(k)));
    nexttile; imshow(pars0.(fname)(:,:,sl) .* mask(:,:,sl), [0 0.2]);
    title(sprintf('S0\\_%d Start', k));
    nexttile; imshow(out.final.(fname)(:,:,sl), [0 0.2]);
    title(sprintf('S0\\_%d Fitted', k));
end
nexttile; imshow(R2s_hat_GT(:,:,sl) .* mask(:,:,sl), [0 50]); title('R2s\_hat GT [1/s]');
nexttile; imshow(pars0.R2s_hat(:,:,sl) .* mask(:,:,sl), [0 50]); title('R2s\_hat Start');
nexttile; imshow(out.final.R2s_hat(:,:,sl), [0 50]); title('R2s\_hat Fitted');

nexttile; imshow(dR2s_dalpha_GT(:,:,sl) .* mask(:,:,sl), [0 10]); title('dR2s\_dalpha GT [1/s/rad]');
nexttile; imshow(pars0.dR2s_dalpha(:,:,sl) .* mask(:,:,sl), [0 10]); title('dR2s\_dalpha Start');
nexttile; imshow(out.final.dR2s_dalpha(:,:,sl), [0 10]); title('dR2s\_dalpha Fitted');
