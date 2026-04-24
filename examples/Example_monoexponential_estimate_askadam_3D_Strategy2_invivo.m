addpath(genpath('../../gacelle/'))
clear

%% Load in vivo MRI data using NIfTI files
% Expected inputs:
%   - one NIfTI file per echo time, e.g. echo1.nii, echo2.nii, ...
%   - one NIfTI mask file, e.g. mask.nii
%
% Edit the paths and echo times below to match your data.

% Echo times in seconds (edit to match your acquisition)
t = [2.5e-3, 5e-3, 10e-3, 15e-3, 20e-3, 25e-3, 30e-3];

% Paths to NIfTI files for each echo (one file per echo time)
echoFiles = { ...
    'path/to/echo1.nii', ...
    'path/to/echo2.nii', ...
    'path/to/echo3.nii', ...
    'path/to/echo4.nii', ...
    'path/to/echo5.nii', ...
    'path/to/echo6.nii', ...
    'path/to/echo7.nii'  ...
};

% Path to the brain mask NIfTI file (values > 0 indicate brain voxels)
maskFile = 'path/to/mask.nii';

% Validate that echo time and file counts match
assert(numel(t) == numel(echoFiles), ...
    'Number of echo times must match number of echo NIfTI files.');

% Load mask
maskVol = niftiread(maskFile);
mask    = maskVol > 0;

% Load echo volumes and stack into a 4D array [Nx, Ny, Nz, Nechoes]
vol1 = niftiread(echoFiles{1});
[Nx, Ny, Nz] = size(vol1);
Nechoes = numel(echoFiles);
y = zeros(Nx, Ny, Nz, Nechoes, 'single');
y(:,:,:,1) = single(vol1);
for k = 2:Nechoes
    y(:,:,:,k) = single(niftiread(echoFiles{k}));
end

%% Set up fitting algorithm
modelParams = {'S0','R2star'};

% Starting point: initialise from a simple ratio of first/last echo for R2*
% and the first echo amplitude for S0
S0init     = y(:,:,:,1);
% guard against log of zero or negative values
ratio      = y(:,:,:,end) ./ max(y(:,:,:,1), eps('single'));
dt         = t(end) - t(1);
R2init     = max(-log(max(ratio, eps('single'))) / dt, 0);

pars0.(modelParams{1}) = double(S0init);
pars0.(modelParams{2}) = double(R2init);

% Fitting options
fitting                     = [];
fitting.modelParams         = modelParams;
fitting.lb                  = [0,   0];                          % lower bound  [S0, R2*]
fitting.ub                  = [max(y(mask>0))*2, 500];           % upper bound  [S0, R2*] (s^-1); must be finite: Inf breaks internal [0,1] normalisation
fitting.iteration           = 4000;
fitting.initialLearnRate    = 0.001;
fitting.lossFunction        = 'l1';
fitting.tol                 = 1e-4;
fitting.convergenceValue    = 1e-8;
fitting.convergenceWindow   = 20;
fitting.isDisplay           = false;
fitting.isOptimiseMemory    = true;

% Forward model (same as the synthetic example)
modelFWD = @Example_monoexponential_FWD_askadam_3D_Strategy2;

% Equal weights across echoes
weights = [];

%% Run optimisation
askadam_obj = askadam;
out = askadam_obj.optimisation(y, mask, weights, pars0, fitting, modelFWD, t, mask);

%% Display results
% Choose a representative slice for 2-D display
sliceIdx = round(Nz / 2);

figure; tiledlayout(1,2);
nexttile; imshow(out.final.S0(:,:,sliceIdx)  .* mask(:,:,sliceIdx), []); title('S0 Fitted');   colorbar;
nexttile; imshow(out.final.R2star(:,:,sliceIdx) .* mask(:,:,sliceIdx), []); title('R2* Fitted (s^{-1})'); colorbar;
