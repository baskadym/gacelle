addpath(genpath('/myriadfs/home/rmapkdy/Scratch/gacelle'))
clear

%% Load in vivo MRI data using NIfTI files
% Expected inputs:
%   - one NIfTI file per echo time, e.g. echo1.nii, echo2.nii, ...
%   - one NIfTI mask file, e.g. mask.nii
%
% Edit the paths and echo times below to match your data.

% Echo times in seconds (edit to match your acquisition)
t = [2.3:2.38:14.2]*1e-3;

% Paths to NIfTI files for each echo (one file per echo time)
echoFiles = dir('/myriadfs/home/rmapkdy/Scratch/input/20220504.M700350/MORSE_v15.8/pdw_scan1/mag/mag*');

% Validate that echo time and file counts match
assert(numel(t) == numel(echoFiles), ...
    'Number of echo times must match number of echo NIfTI files.');

% Load echo volumes and stack into a 4D array [Nx, Ny, Nz, Nechoes]
vol1 = niftiread(fullfile(echoFiles(1).folder, echoFiles(1).name));
[Nx, Ny, Nz] = size(vol1);
Nechoes = numel(echoFiles);
y = zeros(Nx, Ny, Nz, Nechoes, 'single');
y(:,:,:,1) = single(vol1);
for k = 2:Nechoes
    y(:,:,:,k) = single(niftiread(fullfile(echoFiles(k).folder, echoFiles(k).name)));
end
mask = niftiread(fullfile(echoFiles(k).folder, 'mask.nii') );
mask(isinf(mask))=0;
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
fitting.lb                  = [0,   0];    % lower bound  [S0, R2*]
fitting.ub                  = [2*max(S0init(:)), 500];  % upper bound  [S0, R2*] (s^-1)
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
nexttile; imshow(out.final.S0(:,:,sliceIdx)  .* mask(:,:,sliceIdx)); title('S0 Fitted');   colorbar;
nexttile; imshow(out.final.R2star(:,:,sliceIdx) .* mask(:,:,sliceIdx)); title('R2* Fitted (s^{-1})'); colorbar;
save( '/home/rmapkdy/Scratch/output/fitted_para.mat', 'out')
saveas(gcf, '/home/rmapkdy/Scratch/output/gacelle_test_invivo.png')
