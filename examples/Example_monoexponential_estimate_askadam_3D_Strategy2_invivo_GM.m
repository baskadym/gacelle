addpath(genpath('/myriadfs/home/rmapkdy/Scratch/gacelle'))
clear

%% Load two sets of multi-echo in vivo MRI data using NIfTI files
% Expected inputs:
%   - one NIfTI file per echo time for each scan, e.g. mag001.nii, mag002.nii, ...
%   - one NIfTI mask file, e.g. mask.nii
%   - two flip angles, one per scan
%
% The two scans may have different numbers of echoes.
% Edit the paths, echo times and flip angles below to match your data.

% Flip angles in degrees (one per scan)
alpha_values = [6, 24];

% Echo times in seconds for each scan (edit to match your acquisition)
t1 = [2.3:2.38:14.2]*1e-3;   % scan 1 (PDW, alpha = 6 deg)
t2 = [2.3:2.38:14.2]*1e-3;   % scan 2 (T1W, alpha = 24 deg)

% NIfTI file lists for each scan (one file per echo time)
echoFiles1 = dir('/myriadfs/home/rmapkdy/Scratch/input/20220504.M700350/MORSE_v15.8/pdw_scan1/mag/mag*');
echoFiles2 = dir('/myriadfs/home/rmapkdy/Scratch/input/20220504.M700350/MORSE_v15.8/t1w_scan1/mag/mag*');

% Validate that echo time and file counts match for each scan
assert(numel(t1) == numel(echoFiles1), ...
    'Scan 1: number of echo times must match number of echo NIfTI files.');
assert(numel(t2) == numel(echoFiles2), ...
    'Scan 2: number of echo times must match number of echo NIfTI files.');

%% Load echo volumes and stack into 4D arrays [Nx, Ny, Nz, Nechoes]
vol1 = niftiread(fullfile(echoFiles1(1).folder, echoFiles1(1).name));
[Nx, Ny, Nz] = size(vol1);

Nechoes1 = numel(echoFiles1);
y1 = zeros(Nx, Ny, Nz, Nechoes1, 'single');
y1(:,:,:,1) = single(vol1);
for k = 2:Nechoes1
    y1(:,:,:,k) = single(niftiread(fullfile(echoFiles1(k).folder, echoFiles1(k).name)));
end

Nechoes2 = numel(echoFiles2);
y2 = zeros(Nx, Ny, Nz, Nechoes2, 'single');
for k = 1:Nechoes2
    y2(:,:,:,k) = single(niftiread(fullfile(echoFiles2(k).folder, echoFiles2(k).name)));
end

% Load brain mask (shared between the two scans; assumed to be co-registered)
mask = niftiread(fullfile(echoFiles1(1).folder, 'mask.nii'));
mask(isinf(mask)) = 0;

% Concatenate both scans along the echo (4th) dimension: [Nx, Ny, Nz, Nechoes1+Nechoes2]
y = cat(4, y1, y2);

% Build concatenated echo-time and flip-angle vectors (one entry per measurement)
t     = [t1(:)', t2(:)'];
alpha = [alpha_values(1)*ones(1, Nechoes1), alpha_values(2)*ones(1, Nechoes2)];

%% Set up fitting algorithm
modelParams = {'S0', 'R2s_hat', 'dR2s_dalpha'};

% Starting point: initialise S0 and R2s_hat from the PDW scan (scan 1)
S0init  = y1(:,:,:,1);
ratio   = y1(:,:,:,end) ./ max(y1(:,:,:,1), eps('single'));
dt      = t1(end) - t1(1);
R2init  = max(-log(max(ratio, eps('single'))) / dt, 0);

pars0.S0          = double(S0init);
pars0.R2s_hat     = double(R2init);
pars0.dR2s_dalpha = zeros(Nx, Ny, Nz);   % dR2*/dalpha (s^-1 / degree)

% Fitting options
fitting                     = [];
fitting.modelParams         = modelParams;
fitting.lb                  = [0,   0,   -5];   % lower bounds  [S0, R2s_hat, dR2s_dalpha]
fitting.ub                  = [2*max(S0init(:)), 500, 5];  % upper bounds
fitting.iteration           = 4000;
fitting.initialLearnRate    = 0.001;
fitting.lossFunction        = 'l1';
fitting.tol                 = 1e-4;
fitting.convergenceValue    = 1e-8;
fitting.convergenceWindow   = 20;
fitting.isDisplay           = false;
fitting.isOptimiseMemory    = true;

% Forward model for the two-flip-angle monoexponential decay
% S = S0 * exp(-t * (R2s_hat + dR2s_dalpha * alpha))
modelFWD = @Example_monoexponential_FWD_askadam_3D_Strategy2_GM;

% Equal weights across echoes
weights = [];

%% Run optimisation
askadam_obj = askadam;
out = askadam_obj.optimisation(y, mask, weights, pars0, fitting, modelFWD, t, alpha, mask);

%% Display results
% Choose a representative slice for 2-D display
sliceIdx = round(Nz / 2);

figure; tiledlayout(1,3);
nexttile; imshow(out.final.S0(:,:,sliceIdx)          .* mask(:,:,sliceIdx)); title('S0 Fitted');              colorbar;
nexttile; imshow(out.final.R2s_hat(:,:,sliceIdx)     .* mask(:,:,sliceIdx)); title('R2s\_hat Fitted (s^{-1})'); colorbar;
nexttile; imshow(out.final.dR2s_dalpha(:,:,sliceIdx) .* mask(:,:,sliceIdx)); title('dR2s\_dalpha Fitted (s^{-1}/deg)'); colorbar;
save( '/home/rmapkdy/Scratch/output/fitted_para_GM.mat', 'out')
saveas(gcf, '/home/rmapkdy/Scratch/output/gacelle_test_invivo_GM.png')
