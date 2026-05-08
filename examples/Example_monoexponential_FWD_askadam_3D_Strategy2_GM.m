%% S = Example_monoexponential_FWD_askadam_3D_Strategy2_GM(pars, t, alpha, mask)
%
% Input
% --------------
% pars              : input model parameter structure (This is ALWAYS the first input variable)
%   .S0_k           : signal offset at TE=0 for the k-th flip angle (one field per flip angle,
%                     e.g. .S0_1, .S0_2, ... — no parametric form of alpha-dependency assumed)
%   .R2s_hat        : baseline R2* [1/s], common across all flip angles
%   .dR2s_dalpha    : flip-angle sensitivity of R2* [1/s/rad], common across all flip angles
% t                 : [1xNt] echo times [s]
% alpha             : [1xNa] flip angles [rad]
% mask              : 3D signal mask
%
% Output
% --------------
% S                 : signal, [Na*Nt x Nvoxel] (when masked) or [Nx x Ny x Nz x Na*Nt] (otherwise)
%
% Description:
%   Forward model for multi-flip-angle multi-echo monoexponential decay.
%   S0 is flip-angle-dependent: each flip angle k has its own free parameter S0_k,
%   with no assumption on the functional form of S0(alpha).
%   R2s_hat and dR2s_dalpha are shared (common) across all flip angles,
%   giving a common parametric R2* model: R2*(alpha) = R2s_hat + dR2s_dalpha * alpha.
%
%   Signal model:
%       S(alpha_k, TE) = S0_k .* exp(-TE .* (R2s_hat + dR2s_dalpha .* alpha_k))
%
%   Data layout: measurements from all flip angles are concatenated along the
%   4th dimension in the order [alpha_1_TE_1, ..., alpha_1_TE_Nt,
%                                alpha_2_TE_1, ..., alpha_2_TE_Nt, ...].
%
% Kwok-Shing Chan @ MGH
% kchan2@mgh.harvard.edu
% Date created: 7 May 2025
% Date last modified:
%
function S = Example_monoexponential_FWD_askadam_3D_Strategy2_GM(pars, t, alpha, mask)

Na = numel(alpha);

% Put echo times in the 4th dimension for broadcasting
t = reshape(t(:), 1, 1, 1, numel(t));

R2s_hat     = pars.R2s_hat;
dR2s_dalpha = pars.dR2s_dalpha;

% Compute signal for each flip angle and concatenate along the measurement
% dimension (4th dim). No parametric form of S0(alpha) is assumed.
S = [];
for k = 1:Na
    S0_k  = pars.(sprintf('S0_%d', k));
    R2s_k = R2s_hat + dR2s_dalpha .* alpha(k);         % effective R2* at alpha_k
    S_k   = S0_k .* exp(-t .* R2s_k);                  % [Nx,Ny,Nz,Nt] or [1,Nvoxel,1,Nt]
    S     = cat(4, S, S_k);
end

% S is now [Nx,Ny,Nz,Na*Nt] or [1,Nvoxel,1,Na*Nt].
% When parameters are masked (isOptimiseMemory=true) the spatial dims will
% differ from mask; convert to 2D [Na*Nt x Nvoxel] so askAdam can handle it
% without an additional masking step (Strategy 2 memory-efficient approach).
if any(size(R2s_hat, 1:3) ~= size(mask, 1:3))
    S = utils.reshape_ND2GD(S, []);
end

end
