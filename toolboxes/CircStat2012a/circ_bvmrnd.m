function alpha = circ_bvmrnd(params, n)

% alpha = circ_bvmrnd(theta, kappa, n)
%   Simulates n random pairs of angles from a bivariate von Mises distribution,
%   cosine with positive interation with parameters params 
%   
%
%   Input:
%     [params   struct with 5 fields: mu, nu, kappa1, kappa2, kappa3
%     [n        number of samples, default is 10]
%
%
%   Output:
%     alpha     size n x 2: samples from bivariate von Mises distribution
%     cosine with positive interaction
%
%
%   References:
%     Statistics of bivariate von Mises distributions - Marida and Frellsen
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens and Marc J. Velasco, 2009
% velasco@ccs.fau.edu


% default parameter
params_default.mu = 0;
params_default.nu = 0;
params_default.kappa1 = 1;
params_default.kappa2 = 1;
params_default.kappa3 = 0.5;

if nargin < 2
    n = 10;
end
 
if nargin < 1 || isempty(params)
    params = params_default;
end

mu = params(1);
nu = params(2);
kappa1 = params(3);
kappa2 = params(4);
kappa3 = params(5);

phi = zeros(n,1);
psi = zeros(n,1);
% if kappa is small, treat as uniform distribution
if kappa1 < 1e-6 || kappa2 < 1e-6
    phi(:) = 2*pi*rand(n,1);
    psi(:) = 2*pi*rand(n,1);
    alpha = [phi psi];
    return
end

if abs(kappa3) < 1e-6 % no interaction, two independent variables phi and psi
    phi(:) = circ_vmrnd(mu, kappa1, n);
    psi(:) = circ_vmrnd(nu, kappa2, n);
    alpha = [phi psi];
    return
end

% other cases
% simulat psi_p (prime) from the marginal density f(psi) and simulate phi_p from the
% conditional density f(phi|psi = psi_p)

% since we use cosine with positive interaction, the unimodal/bimodal
% conditions follow this model
funA = @(kappa) besseli(1,kappa)./besseli(0,kappa);
fun_kappa13 = @(ang) sqrt(kappa1.^2+kappa3.^2 -2*kappa1.*kappa3.*cos(ang - nu));

%% marginal distribution of psi is unimodal 
if funA(abs(kappa2 - kappa3)) <= abs((kappa2 - kappa3)*kappa1/(kappa2*kappa3))
    % we have to choose the kappa of the von mises such that the distance
    % between marginal density and proposed density is minized. Can I do an
    % fminsearch for this? For now, let's just use kappa2
    psi_star = 0; % 0 offset since the marginal is unimodal
    proposed_marginal_psi = @(ang,kappa) 0.5*circ_vmpdf(ang,nu - psi_star,kappa)+...
    0.5*circ_vmpdf(ang,nu + psi_star,kappa);
    pdf_difference = @(ang,kappa) sum(abs(marginalProb(ang,params)...
        - proposed_marginal_psi(ang,kappa)));
    fminsearch_opts = optimset('PlotFcns',@optimplotfval,'Display','iter',...
        'MaxIter',100,'TolFun',0.01,'TolX',0.01);
    xx = -pi:0.1:pi;
    [kappa_star, ~,~] = fminsearch(@(kappa_star) ...
                pdf_difference(xx,kappa_star),kappa2, fminsearch_opts);
    yy1 = marginalProb(xx,params);
    yy2 = proposed_marginal_psi(xx,kappa_star);
    figure; plot(xx,yy1,'k.',xx,yy2,'r--');
    
    L = threshold_candidate([psi_star kappa_star],params);% threshold L (look at the reference)
    for j = 1:n
        psi_candidate = circ_vmrnd(nu, kappa_star, 1); % candidate
        % threshold L (look at the reference)
        u = rand(1);
        while u > marginalProb(psi_candidate, params)/...
               (L*proposed_marginal_psi(psi_candidate,kappa_star)) % rejection criterion
           psi_candidate = circ_vmrnd(nu, kappa_star, 1);
           u = rand(1);
        end
        psi(j) = psi_candidate;
        psi_nu = atan( - kappa3 *sin(psi_candidate - nu)/(kappa1 - kappa3*cos(psi_candidate-nu)));
        kappa13 = fun_kappa13(psi_candidate);
        phi(j) = circ_vmrnd(psi_nu + mu, kappa13, 1);% generate phi from the conditional distribution phi|psi       
    end
%%  marginal distribution of psi is bimodal 
else
    % theorem 6.4
    fun_psi_star = @(ang) kappa2*kappa3*funA(fun_kappa13(ang))/fun_kappa13(ang) - kappa1;
    syms myangle
    psi_star = vpasolve(fun_psi_star(myangle) == 0,myangle,[0 pi]);
    psi_star = double(psi_star);
    proposed_marginal_psi = @(ang,kappa) 0.5*circ_vmpdf(ang,nu - psi_star,kappa)+...
    0.5*circ_vmpdf(ang,nu + psi_star,kappa);
    pdf_difference = @(ang,kappa) sum(abs(marginalProb(ang,params)...
        - proposed_marginal_psi(ang,kappa)));
    fminsearch_opts = optimset('PlotFcns',@optimplotfval,'Display','iter',...
        'MaxIter',100,'TolFun',0.01,'TolX',0.01);
    xx = -pi:0.1:pi;
    [kappa_star, ~,~] = fminsearch(@(kappa_star) ...
                pdf_difference(xx,kappa_star),kappa2, fminsearch_opts);
    L = threshold_candidate([psi_star kappa_star],params);% threshold L (look at the reference)
    yy1 = marginalProb(xx,params);
    yy2 = proposed_marginal_psi(xx,kappa_star);
    figure; plot(xx,yy1,'k.',xx,yy2,'r--');
    for j = 1:n
        v = rand; % equal mixture model
        if v < 0.5
           psi_candidate = circ_vmrnd(nu-psi_star, kappa_star, 1); % candidate
           %proposed_marginal_psi = @(ang,kappa)circ_vmpdf(ang,nu - psi_star,kappa);
        else
           psi_candidate = circ_vmrnd(nu+psi_star, kappa_star, 1); % candidate
           %proposed_marginal_psi = @(ang,kappa)circ_vmpdf(ang,nu + psi_star,kappa);
        end
        u = rand(1);
        while u > marginalProb(psi_candidate, params)/...
               (L*proposed_marginal_psi(psi_candidate,kappa_star)) % rejection criterion
            v = rand; % equal mixture model
            if v < 0.5
                psi_candidate = circ_vmrnd(nu-psi_star, kappa_star, 1); % candidate
                %proposed_marginal_psi = @(ang,kappa)circ_vmpdf(ang,nu - psi_star,kappa);
            else
                psi_candidate = circ_vmrnd(nu+psi_star, kappa_star, 1); % candidate
                %proposed_marginal_psi = @(ang,kappa)circ_vmpdf(ang,nu+psi_star,kappa);
            end
        end
        % generate phi from the conditional distribution phi|psi
        psi(j) = psi_candidate;
        psi_nu = atan( - kappa3 *sin(psi_candidate - nu)/(kappa1 - kappa3*cos(psi_candidate-nu)));
        kappa13 = fun_kappa13(psi_candidate);
        phi(j) = circ_vmrnd(double(psi_nu) + mu, kappa13, 1);       
    end
end
alpha = [phi psi];
figure; ndhist(phi,psi,'axis',[-pi pi -pi pi],'filter','bins',1,'columns');
xlim([-pi pi]); ylim([-pi pi]); xlabel('phi');ylabel('psi');
set(gcf,'color','white') % White background for the figure.
[xxx,yyy] = meshgrid(-pi:0.1:pi,-pi:0.1:pi);
ppp = circ_bvmpdf(xxx,yyy,mu,nu,kappa1,kappa2,kappa3);
ppp = reshape(ppp,size(xxx));
figure;contour(xxx,yyy,ppp,'ShowText','on');axis square;axis tight;
set(gcf,'color','white'); xlabel('phi');ylabel('psi');
end

function [p] = marginalProb(ang, params_marginal)
    %mu = params_marginal(1); 
    nu = params_marginal(2); 
    kappa1 = params_marginal(3); kappa2 = params_marginal(4); kappa3 = params_marginal(5);
    fun = @(x) 2*pi*besseli(0,sqrt(kappa1.^2+kappa3.^2 ...
    -2*kappa1.*kappa3.*cos(x - nu))).*exp(kappa2.*cos(x-nu));
    Cc = integral((@(x)fun(x)),0,2*pi);
    p = Cc^-1*2*pi*besseli(0,sqrt(kappa1.^2+kappa3.^2 ...
    -2*kappa1.*kappa3.*cos(ang - nu))).*exp(kappa2.*cos(ang - nu));
    %p = Cc^-1*fun(ang);
    p = p';
end

function [L] = threshold_candidate(params_proposed,params_marginal)
    %mu = params_marginal(1); 
    nu = params_marginal(2); 
    %kappa1 = params_marginal(3); kappa2 = params_marginal(4);
    %kappa3 = params_marginal(5);
    psi_star = params_proposed(1); %offset
    kappa_star = params_proposed(2);
    proposed_marginal_psi = @(ang) 0.5*circ_vmpdf(ang,nu - psi_star,kappa_star)+...
    0.5*circ_vmpdf(ang,nu + psi_star,kappa_star);
    x = -pi:0.1:pi;
    y1 = marginalProb(x,params_marginal);
    y2 = proposed_marginal_psi(x);
    ind = y2 > 1e-3; % avoid dividing by 0
    ratio_distributions = y1(ind)./y2(ind);
    L = max(ratio_distributions);
    if L <= 1; L = 1.2; end;
    if L > 5; L = 5; end;
end


