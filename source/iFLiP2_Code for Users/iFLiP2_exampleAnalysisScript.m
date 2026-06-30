% 1. read data
Aout = iFLiP2_readData('testData001.iFLiP2');


% 2. calculate MPET with parameters.
SPCRange = [0.4,12.4]; % usually SPCRange(1) ~ t0-0.6; SPCRange(2) ~ as large as possible [constrained by next pulse time] without the edge artifact.
t0 = 1.1512; % this can come from fitting using the software (or from step 3).
BGLTCurve = [];
afterPulseRatio = 0; 
% afterPulseRatio is the constant background [at ns scale] proportional to 
% the macroscopic intensity. It shows at a fraction. 
% Typically ~0.04 [PDA44 detector setting 3] or ~0.08 [PDA44 setting 4] 
% Note that the example data has no afterPulse (simulated data).

[MPET, correctedData] = iFLiP2_calculateMPET(Aout.data, Aout.header, SPCRange, t0, BGLTCurve, afterPulseRatio);
% note: correctedData is the data with background subtraction.



% 3. fit data
fig_on = 1; % or 0 if you don't want figure.
LTIntensity = sum(correctedData(:,:,1), 2);

% 3a. single exp fit
beta0 = []; % or your preferred initial condition. 1x7 vector: [pop1, tau1, t0, IRFSigma, BG];
isFix = []; % or specifying which property to be fix: e.g. [0 1 1 0 1] will fix tau1, t0 and BG.
fittedData1 = h_fitLifetimeBySingleExp(Aout.LTTime, LTIntensity, 12.5, beta0, isFix, fig_on);

% 3b. double exp fit
beta0 = []; % or your preferred initial condition. 1x7 vector: [pop1, tau1, pop2, tau2, t0, IRFSigma, BG];
isFix = []; % or specifying which property to be fix: e.g. [0 1 0 0 1 0 1] will fix tau1, t0 and BG.
fittedData2 = h_fitLifetimeByDoubleExp(Aout.LTTime, LTIntensity, 12.5, beta0, isFix, fig_on);