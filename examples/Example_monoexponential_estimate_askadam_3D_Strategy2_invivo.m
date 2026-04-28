addpath(genpath('../../gacelle/'))
clear

%% Load in vivo multi-echo data
% Directory containing multi-echo NIfTI files (one .nii/.nii.gz file per echo,
% sorted alphabetically so that the sort order matches echo order)
data_dir = '/myriadfs/home/rmapkdy/Scratch/input/20220504.M700350/MORSE_v15.8/pdw_scan1/mag/';

% Load all NIfTI files from the directory
nii_files = dir(fullfile(data_dir, '*.nii*'));
[~, sort_idx] = sort({nii_files.name});
nii_files = nii_files(sort_idx);

% Stack echoes into a 4-D array: [Nx, Ny, Nz, Necho]
y = [];
for kecho = 1:numel(nii_files)
    echo_data = niftiread(fullfile(nii_files(kecho).folder, nii_files(kecho).name));
    y = cat(4, y, single(echo_data));
end

% Echo times in seconds -- update to match your acquisition
% e.g. t = linspace(2.1e-3, 32.0e-3, numel(nii_files));
t = [];   % <-- REQUIRED: fill in echo times before running this script
assert(~isempty(t), 'Echo times ''t'' must be specified before running this script.');

% Brain mask: keep voxels above 10 % of the maximum signal in the first echo
mask = y(:,:,:,1) > 0.1 * max(y(:,:,:,1), [], 'all');

%% Set up fitting algorithm
modelParams = {'S0','R2star'};

[Nx, Ny, Nz] = deal(size(y,1), size(y,2), size(y,3));

% Starting points
pars0.(modelParams{1}) = ones(Nx, Ny, Nz);          % S0
pars0.(modelParams{2}) = 30 * ones(Nx, Ny, Nz);     % R2* (s^-1)

% Fitting options
fitting                     = [];
fitting.modelParams         = {'S0','R2star'};
fitting.lb                  = [0,   0  ];   % lower bounds
fitting.ub                  = [4, 200  ];   % upper bounds
fitting.iteration           = 4000;
fitting.initialLearnRate    = 0.001;
fitting.lossFunction        = 'l1';
fitting.tol                 = 1e-4;
fitting.convergenceValue    = 1e-8;
fitting.convergenceWindow   = 20;
fitting.isDisplay           = false;
fitting.isOptimiseMemory    = true;

% Forward model handle (monoexponential decay, Strategy 2)
modelFWD = @Example_monoexponential_FWD_askadam_3D_Strategy2;

% Equal weights
weights = [];

%% Run optimisation
askadam_obj = askadam;
out         = askadam_obj.optimisation(y, mask, weights, pars0, fitting, modelFWD, t, mask);

%% Display results
figure; tiledlayout(1,2)
nexttile; imshow(out.final.S0(:,:,round(Nz/2)),     []);  title('S0 Fitted')
nexttile; imshow(out.final.R2star(:,:,round(Nz/2)), []);  title('R2* Fitted (s^{-1})')
