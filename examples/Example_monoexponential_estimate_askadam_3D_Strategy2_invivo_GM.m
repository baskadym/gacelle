addpath(genpath('/myriadfs/home/rmapkdy/Scratch/gacelle'))
clear

%% Load N sets of multi-echo in vivo MRI data using NIfTI files
% Expected inputs:
%   - one NIfTI file per echo time for each scan, e.g. mag001.nii, mag002.nii, ...
%   - one NIfTI mask file in the first scan directory, e.g. mask.nii
%   - one flip angle per scan (any number of scans >= 2 supported)
%
% Each scan may have a different number of echoes.
% Edit the paths, echo times and flip angles below to match your data.

% Flip angles in degrees (one per scan)
alpha_values = [6, 9, 26];

% Echo times in seconds for each scan (one row / cell entry per scan)
% Edit to match your acquisition; scans may have different numbers of echoes.
echo_times{1} = [2.3:2.38:14.2]*1e-3;   % scan 1 (alpha = 6 deg)
echo_times{2} = [2.3:2.38:14.2]*1e-3;   % scan 2 (alpha = 9 deg)
echo_times{3} = [2.3:2.38:14.2]*1e-3;   % scan 3 (alpha = 26 deg)

% Directories containing NIfTI echo files for each scan
scan_dirs{1} = '/myriadfs/home/rmapkdy/Scratch/input/t1w_mfc_3dflash_v3i_6deg_180us_0012';
scan_dirs{2} = '/myriadfs/home/rmapkdy/Scratch/input/t1w_mfc_3dflash_v3i_9deg_180us_0011';
scan_dirs{3} = '/myriadfs/home/rmapkdy/Scratch/input/t1w_mfc_3dflash_v3i_26deg_180us_0007';

%% Validate inputs
Nscans = numel(alpha_values);
assert(numel(echo_times) == Nscans, 'echo_times must have one entry per flip angle.');
assert(numel(scan_dirs)  == Nscans, 'scan_dirs must have one entry per flip angle.');

echoFiles = cell(1, Nscans);
for s = 1:Nscans
    echoFiles{s} = dir(fullfile(scan_dirs{s}, '*.nii'));
    assert(numel(echo_times{s}) == numel(echoFiles{s}), ...
        sprintf('Scan %d: number of echo times must match number of echo NIfTI files.', s));
end

%% Load echo volumes and stack into 4D arrays [Nx, Ny, Nz, Nechoes]
% Determine spatial dimensions from the first echo of the first scan
vol1 = niftiread(fullfile(echoFiles{1}(1).folder, echoFiles{1}(1).name));
[Nx, Ny, Nz] = size(vol1);

y_scans = cell(1, Nscans);
for s = 1:Nscans
    Nechoes_s = numel(echoFiles{s});
    y_s = zeros(Nx, Ny, Nz, Nechoes_s, 'single');
    for k = 1:Nechoes_s
        y_s(:,:,:,k) = single(niftiread(fullfile(echoFiles{s}(k).folder, echoFiles{s}(k).name)));
    end
    y_scans{s} = y_s;
end

% Load brain mask (assumed to be co-registered; located in the first scan directory)
mask = niftiread(fullfile(echoFiles{1}(1).folder, 'mask.nii'));
mask(isinf(mask)) = 0;

% Concatenate all scans along the echo (4th) dimension
y = cat(4, y_scans{:});

% Build concatenated echo-time, flip-angle and scan-index vectors
% (one entry per measurement)
t        = [];
alpha    = [];
scan_idx = [];
for s = 1:Nscans
    t        = [t,        echo_times{s}(:)'];                               %#ok<AGROW>
    alpha    = [alpha,    alpha_values(s)*ones(1, numel(echo_times{s}))];   %#ok<AGROW>
    scan_idx = [scan_idx, s*ones(1, numel(echo_times{s}))];                 %#ok<AGROW>
end

%% Set up fitting algorithm
modelParams = {'S0', 'R2s_hat', 'dR2s_dalpha'};

% Starting point: initialise one S0 map per scan and the shared decay
% parameters from scan 1
y1 = y_scans{1};
S0init  = zeros(Nx, Ny, Nz, Nscans, 'single');
for s = 1:Nscans
    S0init(:,:,:,s) = y_scans{s}(:,:,:,1);
end
ratio   = y1(:,:,:,end) ./ max(y1(:,:,:,1), eps('single'));
dt      = echo_times{1}(end) - echo_times{1}(1);
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

% Forward model for multi-flip-angle monoexponential decay
% S = S0(scan) * exp(-t * (R2s_hat + dR2s_dalpha * alpha))
modelFWD = @Example_monoexponential_FWD_askadam_3D_Strategy2_GM;

% Equal weights across echoes
weights = [];

%% Run optimisation
askadam_obj = askadam;
out = askadam_obj.optimisation(y, mask, weights, pars0, fitting, modelFWD, t, alpha, scan_idx, mask);

%% Display results
% Choose a representative slice for 2-D display
sliceIdx = round(Nz / 2);

figure; tiledlayout(1,Nscans+2);
for s = 1:Nscans
    nexttile;
    imshow(out.final.S0(:,:,sliceIdx,s) .* mask(:,:,sliceIdx));
    title(sprintf('S0 Fitted (\\alpha=%g^\\circ)', alpha_values(s)));
    colorbar;
end
nexttile; imshow(out.final.R2s_hat(:,:,sliceIdx)     .* mask(:,:,sliceIdx)); title('R2s\_hat Fitted (s^{-1})'); colorbar;
nexttile; imshow(out.final.dR2s_dalpha(:,:,sliceIdx) .* mask(:,:,sliceIdx)); title('dR2s\_dalpha Fitted (s^{-1}/deg)'); colorbar;
save( '/home/rmapkdy/Scratch/output/fitted_para_GM.mat', 'out')
saveas(gcf, '/home/rmapkdy/Scratch/output/gacelle_test_invivo_GM.png')
