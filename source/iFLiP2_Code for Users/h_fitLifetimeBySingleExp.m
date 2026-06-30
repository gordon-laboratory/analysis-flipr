function fittedData = h_fitLifetimeBySingleExp(t,intensity, pulseInterval, beta0, isFix, fig_on)

% Input:
%   t: n x 1 vector
%   intensity: cumulative counts at t
%   pulseInterval [in ns] = 1/laserSyncRate * 1e9. Default is 12.5ns
%   beta0: initial condition. 1x5 vector:  [pop1, tau1, t0, IRFSigma, BG];
%   isFix correspond to which property to be fix: e.g. [0 1 0 0 1] will fix tau1 and BG.
%
% Example syntax
%   fittedData = h_fitLifetimeBySingleExp(LTTime,LTIntensity);
%   fittedData = h_fitLifetimeBySingleExp(LTTime,LTIntensity, [], [100 1 1 0.13 0.03])
%   fittedData = h_fitLifetimeBySingleExp(LTTime,LTIntensity, [], [100 1 1 0.13 0.03], [0 0 1 1 1])

% 
t_step = mean(diff(t));

t = t(:); % ensure that it is row vector.
intensity = intensity(:);

if ~exist('beta0', 'var') || isempty(beta0)
    beta0(1) = max(intensity); % population
    [~, peakIdx] = max(intensity);
    estimatedT0 = t(peakIdx) - 0.5;
    MPET = sum(t .* intensity)/sum(intensity) - estimatedT0;
    beta0(2) = MPET; % tau
    beta0(3) = estimatedT0; % t0
    beta0(4) = 0.13; % IRF [sigma of Gaussian]
    beta0(5) = 0; % constant background. by default, don't use background
end

if ~exist('isFix', 'var') || isempty(isFix)
    isFix = [0 0 0 0 1]; % by default, don't use background
end

if ~exist('pulseInterval', 'var') || isempty(pulseInterval)
    pulseInterval = 12.5; %ns, corresponding to 80 MHz
end

if ~exist('fig_on', 'var') || isempty(fig_on)
    fig_on = 1;
end

weight = sqrt(intensity)/sqrt(max(intensity));
% weight(lifetime < 1)=1; % this is copied from ss_firexp2gauss. Now think it is an error, should be below:
weight(intensity < 1) = 1; % not 0 or 0.1! see below.
% lower value is higher weight 

global fittingState; % use fittingState to pass fix opts to the fit function. 

fittingState.isFix = isFix; % this is to avoid variable with identical names.
fittingState.beta0 = beta0; % this is to avoid variable with identical names.

fittingState.pulseI = pulseInterval; % pulse interval in ns for calculating the effect of pre-pulse

[betahat,R,J,converge] = ss_nlinfit(t, intensity, weight, @internal_expGauss, beta0);

% note: the fast fitting requires additional toolbox and is very sensitive
% to initial conditions.
% [betahat,R,J,converge] = h_nlinfit_fast(t, intensity, weight, @internal_expGauss, beta0);


fittedData.originalXdata = t;
fittedData.originalYdata = intensity;
fittedData.originalTime = t; % x is index, while t is the time.
fittedData.fixOpt.tau1 = isFix(2);
fittedData.fixOpt.t0 = isFix(3);
fittedData.fixOpt.IRFSigma = isFix(4);
fittedData.fixOpt.BG = isFix(5);

x2 = t(1):0.1*t_step:t(end);%10X oversample for fittedData
fittingState.isFix = [0 0 0 0 0];%disable fixing here to perform evaluation using fitted betahat.

fittedData.fittedXdata = x2;
fittedData.fittedYdata = internal_expGauss(betahat, x2');
fittedData.fittedTime = x2;
% betahat((2)) = betahat([2, 4, 5, 6])*fittingState.psPerUnit/1000;
fittedData.pop1 = betahat(1);
fittedData.tau1 = betahat(2);%*fittingState.psPerUnit/1000;
fittedData.pop2 = 0;
fittedData.tau2 = 0;
fittedData.t0 = betahat(3);%*fittingState.psPerUnit/1000;
fittedData.IRFSigma = betahat(4);%*fittingState.psPerUnit/1000;
fittedData.BG = betahat(5);%*fittingState.psPerUnit/1000;
fittedData.beta = betahat;

fittedData.method = 'Single Exp';
fittedData.converge = converge;

if fig_on
    figure, plot(fittedData.originalXdata, fittedData.originalYdata, 'b-');
    set(gca, 'YScale', 'log');
    hold on, plot(fittedData.fittedXdata, fittedData.fittedYdata, 'r-');
end

clear global fittingState






%%%%%%%%%%%%%%%%
function y = internal_expGauss(beta, x)
%beta0(1): peak
%beta0(2): exp tau
%beta0(5): center
%beta0(6): gaussian width
% 1/2*erfc[(s^2-tau*x)/{sqrt(2)*s*tau}] * exp[s^2/2/tau^2 - x/tau]

global fittingState

pulseI = fittingState.pulseI;

fixTau1 = fittingState.isFix(2);
fixT0 = fittingState.isFix(3);
fixIRFSigma = fittingState.isFix(4);
fixBG = fittingState.isFix(5);

if fixTau1
    tau1 = fittingState.beta0(2);
else
    tau1 = beta(2);
end

if fixT0
    t0 = fittingState.beta0(3);
else
    t0=beta(3);
end

if (fixIRFSigma)
    IRFSigma = fittingState.beta0(4);
else
    IRFSigma = beta(4);
end

if (fixBG)
    bg = fittingState.beta0(5);
else
    bg = beta(5);
end

% 1/2*erfc[(s^2-tau*x)/{sqrt(2)*s*tau}] * exp[s^2/2/tau^2 - x/tau]
% the code below will consider the current pulse and pulses up to 10x previous tau
n = round(tau1 / pulseI * 10);
if n<2
    n = 2; % at least 2.
end
tOffset = -t0:pulseI:(-t0+n*pulseI); % 5 can be tuned as needed.
y1 = beta(1)*exp(IRFSigma^2/2/tau1^2 - (x+tOffset)/tau1); % matlab's unique syntax allows this without doing repmat. But x has to be rows.
y2 = erfc((IRFSigma^2-tau1*(x+tOffset))/(sqrt(2)*tau1*IRFSigma));

y=sum(y1.*y2, 2)/2+bg;



%%%%%%%%%%%%%%%%%%%%%%
function [beta,r,J,converge] = ss_nlinfit(X,y,weight,model,beta0)
%Modified by Ryohei (adding weight), then by Haining (add converge output)
% [BETA,R,J] = NLINFIT(X,Y,WEIGHT,FUN,BETA0)
% WEIGHT is weight for fitting. For fitting to photo-measurement,
% you should use sqrt(y).
%
%NLINFIT Nonlinear least-squares data fitting by the Gauss-Newton method.
%   NLINFIT(X,Y,FUN,BETA0) estimates the coefficients of a nonlinear
%   function.  Y is a vector.  X is a vector or matrix with the same
%   number of rows as Y.  FUN is a function that accepts two arguments,
%   a coefficient vector and an array of X values, and returns a vector
%   of fitted Y values.  BETA0 is a vector containing initial guesses for
%   the coefficients.
%
%   [BETA,R,J] = NLINFIT(X,Y,FUN,BETA0) returns the fitted coefficients
%   BETA, the residuals R, and the Jacobian J.  You can use these outputs
%   with NLPREDCI to produce error estimates on predictions, and with
%   NLPARCI to produce error estimates on the estimated coefficients.
%
%   Examples
%   --------
%   FUN can be specified using @:
%      nlintool(x, y, @myfun, b0)
%   where MYFUN is a MATLAB function such as:
%      function yhat = myfun(beta, x)
%      b1 = beta(1);
%      b2 = beta(2);
%      yhat = 1 ./ (1 + exp(b1 + b2*x));
%
%   FUN can also be an inline object:
%      fun = inline('1 ./ (1 + exp(b(1) + b(2)*x))', 'b', 'x')
%      nlintool(x, y, fun, b0)
%
%   See also NLPARCI, NLPREDCI, NLINTOOL.

%   B.A. Jones 12-06-94.
%   Copyright 1993-2000 The MathWorks, Inc. 
% $Revision: 2.20 $  $Date: 2000/05/26 18:53:20 $

%if (nargin<4), error('NLINFIT requires four arguments.'); end
if (nargin<5), error('NLINFIT requires five arguments.'); end

if min(size(y)) ~= 1
   error('Requires a vector second input argument.');
end
y = y(:);
weight = weight(:); %Added by RYohei
%weight = 1;
if size(X,1) == 1 % turn a row vector into a column vector.
   X = X(:);
end

wasnan = (isnan(y) | any(isnan(X),2));
if (any(wasnan))
   y(wasnan) = [];
   X(wasnan,:) = [];
end
n = length(y);

p = length(beta0);
beta0 = beta0(:);

J = zeros(n,p);
beta = beta0;
betanew = beta + 1;
maxiter = 100; %default = 100
iter = 0;
betatol = 0.5E-6;  %default = 1e-4
rtol = 0.5E-6; %default = 1e-4
sse = 1;
sseold = sse;
seps = sqrt(eps);
zbeta = zeros(size(beta));
s10 = sqrt(10);
eyep = eye(p);
zerosp = zeros(p,1);
while (norm((betanew-beta)./(beta+seps)) > betatol | abs(sseold-sse)/(sse+seps) > rtol) & iter < maxiter
   if iter > 0 
      beta = betanew;
   end

   iter = iter + 1;
   yfit = feval(model,beta,X);
   r = y - yfit;
   r = r./weight;  %added by Ryohei
   sseold = r'*r;

   for k = 1:p
      delta = zbeta;
      if (beta(k) == 0)
         nb = sqrt(norm(beta));
         delta(k) = seps * (nb + (nb==0));
      else
         delta(k) = seps*beta(k);
      end
      yplus = feval(model,beta+delta,X);
      J(:,k) = (yplus - yfit)/delta(k);
   end

   Jplus = [J;(1.0E-2)*eyep];
   rplus = [r;zerosp];

   % Levenberg-Marquardt type adjustment 
   % Gauss-Newton step -> J\r
   % LM step -> inv(J'*J+constant*eye(p))*J'*r
   step = Jplus\rplus;
   
   betanew = beta + step;
   yfitnew = feval(model,betanew,X);
   rnew = y - yfitnew;
   rnew = rnew./weight; %Added by RYohei
   sse = rnew'*rnew;
   iter1 = 0;
   while sse > sseold & iter1 < 12
      step = step/s10;
      betanew = beta + step;
      yfitnew = feval(model,betanew,X);
      rnew = y - yfitnew;
      rnew = rnew./weight; %Added by Ryohei
      sse = rnew'*rnew; 
      iter1 = iter1 + 1;
   end
end
if iter == maxiter
    beep;
    beep;
   disp('NLINFIT did NOT converge. Returning results from last iteration.');
   converge = 0;
else
    converge = 1;
end
%disp('Iteration number =');
%disp(iter);

