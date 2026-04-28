addpath(genpath('../../gacelle/'))
clear

%% -----------------------------------------------------------------------
%  Load two sets of multi-echo in vivo data acquired at different flip
%  angles (e.g. PDW scan at alpha=6 deg and T1W scan at alpha=24 deg).
%  The two scans may have different numbers of echoes.
% -----------------------------------------------------------------------

% Flip angles (degrees) for the two scans
alpha_values = [6, 24];

% Directories containing multi-echo NIfTI files (one .nii/.nii.gz per echo)
scan_dirs = { ...
    '/myriadfs/home/rmapkdy/Scratch/input/20220504.M700350/MORSE_v15.8/pdw_scan1/mag/', ...
    '/myriadfs/home/rmapkdy/Scratch/input/20220504.M700350/MORSE_v15.8/t1w_scan1/mag/'  };

% Echo times in seconds for each scan -- update to match your acquisition
% e.g. t_scan{1} = linspace(2.1e-3, 32.0e-3, Nt1);
%      t_scan{2} = linspace(2.1e-3, 12.0e-3, Nt2);
t_scan = { [], [] };   % <-- REQUIRED: fill in echo times for each scan before running this script
assert(~isempty(t_scan{1}) && ~isempty(t_scan{2}), 'Echo times ''t_scan'' must be specified for both scans before running this script.');

%% Load NIfTI files and stack echoes
y_scan = cell(1, numel(scan_dirs));
Nt     = zeros(1, numel(scan_dirs));

for kscan = 1:numel(scan_dirs)
    nii_files = dir(fullfile(scan_dirs{kscan}, '*.nii*'));
    [~, sort_idx] = sort({nii_files.name});
    nii_files = nii_files(sort_idx);
    Nt(kscan) = numel(nii_files);

    tmp = [];
    for kecho = 1:Nt(kscan)
        echo_data = niftiread(fullfile(nii_files(kecho).folder, nii_files(kecho).name));
        tmp = cat(4, tmp, single(echo_data));
    end
    y_scan{kscan} = tmp;
end

% Concatenate both scans along the echo (4th) dimension: [Nx, Ny, Nz, Nt1+Nt2]
y = cat(4, y_scan{1}, y_scan{2});

% Build concatenated echo-time and flip-angle vectors (one entry per measurement)
t     = [t_scan{1}(:)', t_scan{2}(:)'];
alpha = [alpha_values(1) * ones(1, Nt(1)), alpha_values(2) * ones(1, Nt(2))];

[Nx, Ny, Nz] = deal(size(y,1), size(y,2), size(y,3));

%% Brain mask (threshold first echo of first scan)
mask = y_scan{1}(:,:,:,1) > 0.1 * max(y_scan{1}(:,:,:,1), [], 'all');

% Optional: apply an additional gray-matter (GM) tissue mask loaded externally
% gm_mask = niftiread('path/to/gm_mask.nii.gz') > 0;
% mask    = mask & gm_mask;

%% Set up fitting algorithm
modelParams = {'S0', 'R2s_hat', 'dR2s_dalpha'};

% Starting points
pars0.S0          = ones(Nx, Ny, Nz);           % proton-density-weighted amplitude
pars0.R2s_hat     = 30  * ones(Nx, Ny, Nz);     % R2* at alpha=0 (s^-1)
pars0.dR2s_dalpha = zeros(Nx, Ny, Nz);          % dR2*/dalpha (s^-1 / degree)

% Fitting options
fitting                     = [];
fitting.modelParams         = modelParams;
fitting.lb                  = [0,   0,   -5 ];   % lower bounds  [S0, R2s_hat, dR2s_dalpha]
fitting.ub                  = [4, 200,    5 ];   % upper bounds
fitting.iteration           = 4000;
fitting.initialLearnRate    = 0.001;
fitting.lossFunction        = 'l1';
fitting.tol                 = 1e-4;
fitting.convergenceValue    = 1e-8;
fitting.convergenceWindow   = 20;
fitting.isDisplay           = false;
fitting.isOptimiseMemory    = true;

% Forward model handle (two-flip-angle monoexponential, Strategy 2)
modelFWD = @Example_monoexponential_FWD_askadam_3D_Strategy2_GM;

% Equal weights
weights = [];

%% Run optimisation
% Extra inputs after modelFWD are passed to the forward model as:
%   modelFWD(pars, t, alpha, mask)
askadam_obj = askadam;
out         = askadam_obj.optimisation(y, mask, weights, pars0, fitting, modelFWD, t, alpha, mask);

%% Display results
sl = round(Nz/2);   % central axial slice
figure; tiledlayout(1,3)
nexttile; imshow(out.final.S0(:,:,sl),          []);  title('S0 Fitted')
nexttile; imshow(out.final.R2s_hat(:,:,sl),     []);  title('R2s\_hat Fitted (s^{-1})')
nexttile; imshow(out.final.dR2s_dalpha(:,:,sl), []);  title('dR2s\_dalpha Fitted (s^{-1}/deg)')
