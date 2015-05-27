function [p, phi, psi] = circ_bvmpdf(phi,psi, mu, nu, kappa1,  kappa2, kappa3, opts)

% [p, phi, psi] = circ_bvmpdf(phi,psi, mu, nu, kappa1,  kappa2, kappa3)
%   Computes the bivariate circular von Mises pdf with preferred directions
%   mu and nu and concentrations kappa1, kappa2, correlation kappa3
%   at each of the angles in phi and psi 
%
%   The vmpdf is given by f(phi) =
%   Cc*exp{kappa1*cos(phi-mu) + kappa2*cose(psi-nu) - kappa3*cos(phi-mu -
%   psi+nu)}
%
%   Input:
%     phi,psi   angles to evaluate pdf at, if empty phi, psi are chosen to
%               100 uniformly spaced points around the circle
%     mu, nu    preferred direction, default is 0
%     kappa1, kapp2    concentration parameter, default is 1
%     kapp3     correlation, default is (kappa1 + kappa2)/2
%     opts      sine or poscosine model
%   Output:
%     p         von Mises pdf evaluated at phi, psi
%     phi,psi     angles at which pdf was evaluated
%
%
%   References:
%     Statistical analysis of circular data, Fisher
%
% Circular Statistics Toolbox for Matlab

% 1D By Philipp Berens and Marc J. Velasco, 2009
% velasco@ccs.fau.edu

% modified by Luong Nguyen, 2015
% lun5@pitt.edu
% if no angles are supplied, 100 evenly spaced points around the circle are
% chosen
if nargin < 2 || isempty(phi) || isempty(psi)
    phi = linspace(-pi,pi, 101)';
    phi = phi(1:end-1);
    psi = linspace(-pi,pi, 101)';
    psi = psi(1:end-1);
end

if nargin < 3
    mu = 0;
end

if nargin < 4
    nu = 0;
end

if nargin < 5
    kappa1 = 1;
end

if nargin < 6
    kappa2 = 1;
end

if nargin < 7
    kappa3 = (kappa1 + kappa2)/2;
end

if nargin < 8 
    opts.model = 'pcosine';
end

model = opts.model;

%phi = phi(:);
%psi = psi(:);

if kappa1 < 0 || kappa2 < 0
    error('concentration parameters have to be positive');
end

max_kappa = 100; %before it's 150
if kappa1 > max_kappa
    kappa1 = max_kappa + rand;
end

if kappa2 > max_kappa
    kappa2 = max_kappa + rand;
end

% evaluate pdf
if abs(kappa3) < 1e-5 % phi psi independent case
    C = 1/(2*pi*besseli(0,kappa1) * 2*pi*besseli(0,kappa2));
    p = C * exp(kappa1.*cos(phi-mu) + kappa2.*cos(psi-nu) - kappa3.*cos(phi-mu -psi+nu));
    return;
end

if strcmp(model,'pcosine') % positive cosine    
    fun = @(x, nu, kappa1, kappa2, kappa3) 2*pi*besseli(0,sqrt(kappa1.^2+kappa3.^2 ...
    -2*kappa1.*kappa3.*cos(x - nu))).*exp(kappa2.*cos(x-nu));
elseif strcmp(model,'sine') % sine model
    fun = @(x, nu, kappa1, kappa2, kappa3) 2*pi*besseli(0,sqrt(kappa1.^2 + ...
        kappa3^2.*(sin(x - nu)).^2)).* exp(kappa2.*cos(x-nu));
end

C_inv = integral((@(x)fun(x, nu, kappa1, kappa2, kappa3)),0,2*pi);
C = C_inv.^(-1);
p = C .* exp(kappa1.*cos(phi-mu) + kappa2.*cos(psi-nu) - kappa3.*cos(phi-mu -psi+nu));
end
