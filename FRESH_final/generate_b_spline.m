function [phi, t] = generate_b_spline(N, T_s, T, mode)
% -------------------------------------------------------------------------
% Communications and Signal Processing Group
% Department of Electrical and Electronic Engineering
% Imperial College London, 2011
%
% Date        : 21/11/2011
% Supervisor  : Dr Pier Luigi Dragotti
% Author      : Jon Onativia
%
% File        : generate_b_spline.m
% -------------------------------------------------------------------------
% Generate the B-spline of order N.
%   phi(t) = beta_0(t) * beta_0(t) * ... * beta_0(t)
%
% USAGE:
%  [phi, t] = generate_b_spline(N, T_s[, T, mode])
%
% INPUT:
%  - N      : Order of the B-spline.
%  - T_s    : Time resolution of the spline.
%  - T      : Optional argument. Scale factor. Default T = 1.
%  - mode   : Optional argument. 'causal', 'symmetric' or 'anticausal'. 
%             Default mode = 'causal'.
%
% OUTPUT:
%  - phi    : B-Spline.
%  - t      : Time stamps of the corresponding values of the B-spline.
%

if nargin < 2 || nargin > 4
    error('generate_b_spline:err_arg', 'The number of input arguments is incorrect.')
elseif nargin < 4
    mode = 'causal';
    if nargin < 3
        T = 1;
    end
end

% Apply scaling factor
T_s = T_s / T;

t_phi_1          = 0;
t_phi_2          = 1;
t_phi            = (t_phi_1:T_s:t_phi_2)';
sub_phi          = ones(length(t_phi), N+1);

sub_phi(1,:) = 1/2;
sub_phi(end,:) = 1/2;

phi   = sub_phi(:,1);
t_0   = t_phi(1);
t_end = t_phi(end);
for i = 1:N
    t_0   = t_0 + t_phi(1);
    t_end = t_end + t_phi(end);
    phi   = T_s * conv(phi, sub_phi(:,i+1));
end

t = (t_0:T_s:t_end)';
t = t * T;

if strcmp(mode, 'symmetric')
    t_mid      = (t(end) - t(1)) / 2;
    t          = t - t_mid;
    [aaa, i_max] = max(phi);
    if phi(i_max) ~= phi(t == 0)
        t = t - t(i_max);
    end
elseif strcmp(mode, 'anticausal')
    phi = phi(end:-1:1);
    t   = -t(end:-1:1);
end