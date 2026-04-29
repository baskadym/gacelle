%% Example: Joint R1-R2* mapping from in-vivo variable flip angle (vFA) multi-echo GRE data
%
% This script demonstrates how to use gpuJointR1R2starMapping to estimate
% R1, M0 and R2* from in-vivo multi-echo GRE data acquired at an arbitrary
% number of flip angles.  Each flip angle dataset lives in its own folder
% and contains one NIfTI file per echo (or a single 4-D NIfTI with all
% echoes concatenated along the 4th dimension).
%
% ---------------------------------------------------------------------------
% QUICK-START: edit the five blocks marked with  <<< USER INPUT >>>  below.
% ---------------------------------------------------------------------------
%
% Kwok-Shing Chan @ MGH / expanded by GitHub Copilot
% Date created: 2025
% -------------------------------------------------------------------------

addpath(genpath('../../gacelle/'))

%% =========================================================================
%  <<< USER INPUT 1 >>>  Sequence parameters
% ==========================================================================

% Flip angles (degrees) – one value per folder listed in data_dirs
alpha_values = [6, 9, 26];

% Echo times (seconds) – must be identical for every flip-angle acquisition
te = [2.2e-3, 4.6e-3, 7.0e-3, 9.4e-3, 11.8e-3];   % example: 5 echoes

% Repetition time (seconds)
tr = 20e-3;

%% =========================================================================
%  <<< USER INPUT 2 >>>  Data directories
%  One directory per flip angle, in the same order as alpha_values.
%  Each directory must contain the magnitude NIfTI file(s) for that FA.
% ==========================================================================

data_dirs = { ...
    '/myriadfs/home/rmapkdy/Scratch/input/t1w_mfc_3dflash_v3i_6deg_180us_0012', ...
    '/myriadfs/home/rmapkdy/Scratch/input/t1w_mfc_3dflash_v3i_9deg_180us_0011', ...
    '/myriadfs/home/rmapkdy/Scratch/input/t1w_mfc_3dflash_v3i_26deg_180us_0007', ...
};

%% =========================================================================
%  <<< USER INPUT 3 >>>  Brain mask (optional)
%  Set mask_file = '' to auto-generate a mask from the data.
% ==========================================================================

mask_file = '';   % e.g. '/path/to/brain_mask.nii'

%% =========================================================================
%  <<< USER INPUT 4 >>>  B1+ map (optional)
%  A relative B1 map (actual FA / nominal FA).  Set b1_file = '' to assume
%  a spatially uniform B1 (equivalent to b1 = 1 everywhere).
% ==========================================================================

b1_file = '';     % e.g. '/path/to/b1map.nii'

%% =========================================================================
%  <<< USER INPUT 5 >>>  Output filename (optional)
%  Leave empty to skip saving to disk.
% ==========================================================================

output_filename = '';   % e.g. 'results_joint_R1R2star.mat'

%% =========================================================================
%  Load multi-echo magnitude data
% ==========================================================================

% Validate that alpha_values and data_dirs have the same length
if numel(alpha_values) ~= numel(data_dirs)
    error('alpha_values and data_dirs must have the same number of elements.');
end

NFA = numel(alpha_values);
img = [];   % will become [Nx, Ny, Nz, Necho, NFA]

for kfa = 1:NFA

    fprintf('Loading data for flip angle %d deg from:\n  %s\n', alpha_values(kfa), data_dirs{kfa});

    % Collect all NIfTI files in the directory (sorted alphabetically so
    % that echo order is preserved when files are named consistently).
    nii_files = dir(fullfile(data_dirs{kfa}, '*.nii'));
    if isempty(nii_files)
        nii_files = dir(fullfile(data_dirs{kfa}, '*.nii.gz'));
    end
    if isempty(nii_files)
        error('No NIfTI files found in directory:\n  %s', data_dirs{kfa});
    end

    % Sort by filename to ensure correct echo ordering
    [~, sort_idx] = sort({nii_files.name});
    nii_files     = nii_files(sort_idx);

    % Load and concatenate echoes along the 4th dimension
    fa_img = [];
    for kfile = 1:numel(nii_files)
        nii_path = fullfile(nii_files(kfile).folder, nii_files(kfile).name);
        vol      = single(niftiread(nii_path));
        % Support both 3-D (one echo per file) and 4-D (all echoes in one file)
        fa_img   = cat(4, fa_img, vol);
    end

    % Concatenate this flip angle along the 5th dimension
    img = cat(5, img, fa_img);

end

fprintf('Data loaded: [%s] (Nx, Ny, Nz, Necho, NFA)\n', num2str(size(img)));

% Basic consistency check: number of echoes must match te
if size(img, 4) ~= numel(te)
    error(['Number of echoes in loaded data (%d) does not match the number ' ...
           'of echo times specified (%d). Please check ''te'' and your data directories.'], ...
           size(img, 4), numel(te));
end

%% =========================================================================
%  Load (or generate) brain mask
% ==========================================================================

if ~isempty(mask_file)
    mask = logical(niftiread(mask_file));
    fprintf('Brain mask loaded from: %s\n', mask_file);
else
    % Auto-generate mask: voxels whose maximum signal across all echoes and
    % flip angles exceeds 5 % of the global maximum are considered brain.
    mask = max(max(abs(img), [], 4), [], 5) ./ max(abs(img(:))) > 0.05;
    fprintf('Brain mask auto-generated from data (threshold = 5%% of max signal).\n');
end

%% =========================================================================
%  Load B1+ map (optional)
% ==========================================================================

extraData = [];
if ~isempty(b1_file)
    b1           = single(niftiread(b1_file));
    extraData.b1 = b1;
    fprintf('B1 map loaded from: %s\n', b1_file);
else
    fprintf('No B1 map provided; assuming uniform flip angle (b1 = 1).\n');
    % gpuJointR1R2starMapping handles extraData = [] gracefully by setting b1 = 1
end

%% =========================================================================
%  Set up the fitting algorithm
% ==========================================================================

fitting                     = [];
fitting.iteration           = 4000;
fitting.initialLearnRate    = 0.001;
fitting.decayRate           = 0;
fitting.convergenceValue    = 1e-8;
fitting.convergenceWindow   = 20;
fitting.tol                 = 1e-3;
fitting.lossFunction        = 'l1';
fitting.isWeighted          = false;
fitting.outputFilename      = output_filename;

%% =========================================================================
%  Run joint R1–R2* estimation
% ==========================================================================

fprintf('\nStarting joint R1-R2* estimation with %d flip angle(s): [%s] deg\n\n', ...
        NFA, num2str(alpha_values));

obj   = gpuJointR1R2starMapping(te, tr, alpha_values);
out   = obj.estimate(img, mask, extraData, fitting);

%% =========================================================================
%  Display results (central slice)
% ==========================================================================

centre_slice = round(size(mask, 3) / 2);

figure('Name', 'Joint R1-R2* Mapping Results');
tiledlayout(1, 3);

nexttile;
imshow(out.final.M0(:,:,centre_slice) .* mask(:,:,centre_slice), []);
title('M0 (a.u.)');

nexttile;
imshow(out.final.R1(:,:,centre_slice) .* mask(:,:,centre_slice), [0 2]);
title('R1 (s^{-1})');

nexttile;
imshow(out.final.R2star(:,:,centre_slice) .* mask(:,:,centre_slice), [0 100]);
title('R2* (s^{-1})');

colormap gray;
