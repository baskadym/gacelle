%% S = Example_monoexponential_FWD_askadam_3D_Strategy2_GM( pars, t, alpha, scan_idx, mask)
%
% Input
% --------------
% pars          : input model parameter structure (ALWAYS the first input variable)
%   .S0             - proton-density-weighted signal amplitude, one map per scan
%   .R2s_hat        - R2* at flip angle = 0  (s^-1)
%   .dR2s_dalpha    - linear R2* change per unit of flip angle (s^-1 / degree)
% t             : [1 x Nmeas] echo times in seconds, concatenated over all scans
% alpha         : [1 x Nmeas] flip angle (degrees) associated with each measurement
%                 (same length as t, e.g. [alpha1*ones(1,Nt1), alpha2*ones(1,Nt2)])
% scan_idx      : [1 x Nmeas] scan index associated with each measurement
% mask          : 3D signal mask
%
% Output
% --------------
% S             : model signal, [Nmeas x Nvoxel] matrix
%
% Description:
%   Forward model for multi-flip-angle monoexponential decay (Strategy 2).
%   The signal model is:
%       S = S0(scan) .* exp( -t .* (R2s_hat + dR2s_dalpha .* alpha) )
%
% Kwok-Shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 28 April 2025
% Date last modified:
%
function S = Example_monoexponential_FWD_askadam_3D_Strategy2_GM( pars, t, alpha, scan_idx, mask)

% Reshape t, alpha and scan_idx into the 4th dimension for broadcasting with
% N-D spatial arrays (Nx x Ny x Nz x Nmeas)
t        = reshape(t(:),        1, 1, 1, numel(t));
alpha    = reshape(alpha(:),    1, 1, 1, numel(alpha));
scan_idx = reshape(scan_idx(:), 1, 1, 1, numel(scan_idx));

% Retrieve model parameters (N-D spatial arrays)
S0          = pars.S0;
R2s_hat     = pars.R2s_hat;
dR2s_dalpha = pars.dR2s_dalpha;

Nscans = max(scan_idx(:));
assert(size(S0,4) == Nscans, 'pars.S0 must contain one map per scan.');

% Forward signal with scan-specific S0 offsets and common decay parameters
S = zeros([size(S0,1), size(S0,2), size(S0,3), numel(t)], 'like', S0);
for s = 1:Nscans
    idx = reshape(scan_idx == s, 1, []);
    S0_s = S0(:,:,:,s);
    S(:,:,:,idx) = S0_s .* exp( -t(:,:,:,idx) .* (R2s_hat + dR2s_dalpha .* alpha(:,:,:,idx)) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Strategy 2 memory efficiency:
% Inside the askAdam optimisation loop the spatial parameters are masked
% to [1 x Nvoxel], so S has size [1 x Nvoxel x 1 x Nmeas].
% We must convert it to a 2D [Nmeas x Nvoxel] matrix for askAdam.
if any(size(S0,1:3) ~= size(mask,1:3))   % input was masked
    S = utils.reshape_ND2GD(S,[]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
