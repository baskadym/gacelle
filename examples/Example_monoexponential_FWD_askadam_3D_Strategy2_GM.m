%% S = Example_monoexponential_FWD_askadam_3D_Strategy2_GM( pars, t, alpha, mask)
%
% Input
% --------------
% pars          : input model parameter structure (ALWAYS the first input variable)
%   .S0             - proton-density-weighted signal amplitude
%   .R2s_hat        - R2* at flip angle = 0  (s^-1)
%   .dR2s_dalpha    - linear R2* change per unit of flip angle (s^-1 / degree)
% t             : [1 x Nmeas] echo times in seconds, concatenated over all scans
% alpha         : [1 x Nmeas] flip angle (degrees) associated with each measurement
%                 (same length as t, e.g. [alpha1*ones(1,Nt1), alpha2*ones(1,Nt2)])
% mask          : 3D signal mask
%
% Output
% --------------
% S             : model signal, [Nmeas x Nvoxel] matrix
%
% Description:
%   Forward model for two-flip-angle monoexponential decay (Strategy 2).
%   The signal model is:
%       S = S0 .* exp( -t .* (R2s_hat + dR2s_dalpha .* alpha) )
%
% Barbara Dymerska @UCL FIL
%
function S = Example_monoexponential_FWD_askadam_3D_Strategy2_GM( pars, t, alpha, mask)

% Reshape t and alpha into the 4th dimension for broadcasting with
% N-D spatial arrays (Nx x Ny x Nz x Nmeas)
t     = reshape(t(:),     1, 1, 1, numel(t));
alpha = reshape(alpha(:), 1, 1, 1, numel(alpha));

% Retrieve model parameters (N-D spatial arrays)
S0          = pars.S0;
R2s_hat     = pars.R2s_hat;
dR2s_dalpha = pars.dR2s_dalpha;

% Forward signal:  S = S0 * exp( -t * (R2s_hat + dR2s_dalpha * alpha) )
S = S0 .* exp( -t .* (R2s_hat + dR2s_dalpha .* alpha) );

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
