function fittedData = h_fitLifetimeByDoubleExp(t,intensity, pulseInterval, beta0, isFix, fig_on)

% Input:
%   t: n x 1 vector
%   intensity: cumulative counts at t
%   pulseInterval [in ns] = 1/laserSyncRate * 1e9. Default is 12.5ns
%   beta0: initial condition. 1x7 vector: [pop1, tau1, pop2, tau2, t0, IRFSigma, BG];
%   isFix correspond to which property to be fix: e.g. [0 1 0 0 1 0 1] will fix tau1, t0 and BG.
%       Note that pop1 and pop2 will not be fixed (the value has no effect).
%
% Example syntax
%   fittedData = h_fitLifetimeBySingleExp(LTTime,LTIntensity);
%   fittedData = h_fitLifetimeBySingleExp(LTTime,LTIntensity, [], [100 3 200 0.5 1 0.13 0.03])
%   fittedData = h_fitLifetimeBySingleExp(LTTime,LTIntensity, 12.5, [100 3 200 0.5 1 0.13 0.03], [0 1 0 1 1 1 1])


t_step = mean(diff(t));
% fit_range(1) = round(t(1)/t_step);
% fit_range(2) = round(t(end)/t_step);
% x = fit_range(1):fit_range(2); % fitting uses index rather than real time. Don't want to change now but keep in mind.

t = t(:); % ensure that it is row vector.
intensity = intensity(:);

if ~exist('beta0', 'var') || isempty(beta0)
    [maxIntensity, peakIdx] = max(intensity);
    estimatedT0 = t(peakIdx) - 0.5;
    MPET = sum(t .* intensity)/sum(intensity) - estimatedT0;
    beta0(1) = maxIntensity/2;
    beta0(2) = MPET*2;
    beta0(3) = maxIntensity/2;
    beta0(4) = MPET/2;
    beta0(5) = estimatedT0;
    beta0(6) = 0.13;
    beta0(7) = 0;
end

if ~exist('isFix', 'var') || isempty(isFix)
    isFix = [0 0 0 0 0 0 1];
end

if ~exist('fig_on', 'var') || isempty(fig_on)
    fig_on = 1;
end

if ~exist('pulseInterval', 'var') || isempty(pulseInterval)
    pulseInterval = 12.5; %ns, corresponding to 80 MHz
end

weight = sqrt(intensity)/sqrt(max(intensity));
% weight(lifetime < 1)=1; % this is copied from ss_firexp2gauss. Now think it is an error, should be below:
weight(intensity < 1) = 1; % not 0 or 0.1! see below.
% lower value is higher weight 

global fittingState; % use fittingState to pass fix opts to the fit function. 

fittingState.isFix = isFix; % this is to avoid variable with identical names.
fittingState.beta0 = beta0; % this is to avoid variable with identical names.

fittingState.pulseI = pulseInterval; % pulse interval in ns for calculating the effect of pre-pulse

[betahat,R,J,converge] = ss_nlinfit(t, intensity, weight, @internal_exp2Gauss, beta0);


fittedData.originalXdata = t;
fittedData.originalYdata = intensity;
fittedData.originalTime = t; % x is index, while t is the time.
fittedData.fixOpt.tau1 = isFix(2);
fittedData.fixOpt.tau2 = isFix(4);
fittedData.fixOpt.t0 = isFix(5);
fittedData.fixOpt.IRFSigma = isFix(6);
fittedData.fixOpt.BG = isFix(7);

x2 = t(1):0.1*t_step:t(end);%10X oversample for fittedData
fittingState.isFix = [0 0 0 0 0 0 0];%disable fixing here to perform evaluation using fitted betahat.

fittedData.fittedXdata = x2;
fittedData.fittedYdata = internal_exp2Gauss(betahat, x2');
fittedData.fittedTime = x2;
fittedData.pop1 = betahat(1);
fittedData.tau1 = betahat(2);%*fittingState.psPerUnit/1000;
fittedData.pop2 = betahat(3);
fittedData.tau2 = betahat(4);
fittedData.t0 = betahat(5);
fittedData.IRFSigma = betahat(6);
fittedData.BG = betahat(7);
fittedData.beta = betahat;

fittedData.method = 'Double Exp';
fittedData.converge = converge;

if fig_on
    figure, plot(fittedData.originalXdata, fittedData.originalYdata, 'b-');
    set(gca, 'YScale', 'log');
    hold on, plot(fittedData.fittedXdata, fittedData.fittedYdata, 'r-');
end

clear global fittingState





%%%%%%%%%%%%%%%%
function y = internal_exp2Gauss(beta, x)
%beta0(1): peak
%beta0(2): exp tau
%beta0(5): center
%beta0(6): gaussian width
% 1/2*erfc[(s^2-tau*x)/{sqrt(2)*s*tau}] * exp[s^2/2/tau^2 - x/tau]

global fittingState

pulseI = fittingState.pulseI;

fixTau1 = fittingState.isFix(2);
fixTau2 = fittingState.isFix(4);
fixT0 = fittingState.isFix(5);
fixBeta6 = fittingState.isFix(6);
fixBG = fittingState.isFix(7);

if (fixTau1)
    tau1 = fittingState.beta0(2);
else
    tau1 = beta(2);
end

if (fixTau2)
    tau2 = fittingState.beta0(4);
else
    tau2 = beta(4);
end

if fixT0
    t0 = fittingState.beta0(5);
else
    t0=beta(5);
end

if (fixBeta6)
    beta6 = fittingState.beta0(6);
else
    beta6 = beta(6);
end

if (fixBG)
    bg = fittingState.beta0(7);
else
    bg = beta(7);
end

% y1 = beta(1)*exp(beta6^2/2/tau1^2 - (x-t0)/tau1);
% y2 = erfc((beta6^2-tau1*(x-t0))/(sqrt(2)*tau1*beta6));
% y=y1.*y2;
% 
% %"Pre" pulse
% y1 = beta(1)*exp(beta6^2/2/tau1^2 - (x-t0+pulseI)/tau1);
% y2 = erfc((beta6^2-tau1*(x-t0+pulseI))/(sqrt(2)*tau1*beta6));
% 
% ya = y1.*y2+y;
% ya = ya/2;
% 
% y1 = beta(3)*exp(beta6^2/2/tau2^2 - (x-t0)/tau2);
% y2 = erfc((beta6^2-tau2*(x-t0))/(sqrt(2)*tau2*beta6));
% y=y1.*y2;
% 
% y1 = beta(3)*exp(beta6^2/2/tau2^2 - (x-t0+pulseI)/tau2);
% y2 = erfc((beta6^2-tau2*(x-t0+pulseI))/(sqrt(2)*tau2*beta6));
% 
% yb = y1.*y2+y;
% yb = yb/2;
% 
% y=ya+yb;

% basic formula: 1/2*erfc[(s^2-tau*x)/{sqrt(2)*s*tau}] * exp[s^2/2/tau^2 - x/tau]
% 
% above is still limited. Ryohei only do pre pulse. But even after 3x, the
% fitting can benefit from more pre pulses if tau is slow.

n = round(max(abs([tau1, tau2])) / pulseI * 10);
if n<2
    n = 2; % at least 2. Note that 2 is to have pre and pre-pre pulse also considered.
    % when n = 1, the result is equivalent to above.
end
% n = 10;

tOffset = -t0:pulseI:(-t0+n*pulseI); 

y1 = beta(1)*exp(beta6^2/2/tau1^2 - (x+tOffset)/tau1); % matlab's unique syntax allows this without doing repmat. But x has to be rows.
y2 = erfc((beta6^2-tau1*(x+tOffset))/(sqrt(2)*tau1*beta6));

y1a = beta(3)*exp(beta6^2/2/tau2^2 - (x+tOffset)/tau2); % matlab's unique syntax allows this without doing repmat. But x has to be rows.
y2a = erfc((beta6^2-tau2*(x+tOffset))/(sqrt(2)*tau2*beta6));

y=sum(y1.*y2, 2)/2 + sum(y1a.*y2a, 2)/2 + bg;