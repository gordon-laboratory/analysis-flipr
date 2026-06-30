function [MPET, correctedData] = iFLiP2_calculateMPET(data, header, SPCRange, t0, BGLTCurve, afterPulseRatio)

% Calculate MPET from lifetime curves.
%
% INPUTS:
%   data      : lifetime curves, size [nBins x nTime]
%   header    : header structure that contains metadata of acquisition
%   SPCrange  : [tMin tMax], range of LTTime to include, e.g. [0.3 12.3]
%   t0        : in ns, scalar time offset to subtract from MPET
%   BGLTCurve : [Optional] background lifetime curve (single timepoint)
%   afterPulseRatio: [Optional] after pulse ratio. This can be measured, or
%               estimated using the measurements of standard solutions. 
%
% OUTPUT:
%   MPET      : [nTime x nChannels]
%
% Formula:
%   MPET = sum(I(t)*t)/sum(I(t)) - t0
%
% Notes:
%   - Only bins with LTTime >= SPCrange(1) and LTTime <= SPCrange(2) are used.
%   - If total counts = 0, MPET is returned as NaN for that point.

    siz = size(data);
    
    if ~exist('BGLTCurve', 'var') || isempty(BGLTCurve)
        BGLTCurve = zeros(siz(1),1);
    end
    
    samplingFreq = header.state.samplingFreq.Value;
    if isfield(header.init, 'deadTime_ns') % some old data does not have
        deadTime = header.init.deadTime_ns / 1e9;
    else
        deadTime = 25 / 1e9; % 25 ns; assuming TH260p, which is the only device supported by iFLiP2.
    end
    [~, ~, ~, correctedData] = internal_calcEffectiveBG(data, BGLTCurve, deadTime, samplingFreq, afterPulseRatio);
    LTTime = (0:header.acq.LTCurveLength-1)'*header.acq.LTResolution;
    % choose bins within SPC range
    idx = (LTTime >= SPCRange(1)) & (LTTime <= SPCRange(2));
    
    if ~any(idx)
        error('No LTTime bins fall within SPCrange.');
    end
    
    LTTime_use = LTTime(idx);
    data_use = correctedData(idx, :);
    
    MPET = nan(siz(2), 1);
    
    totalCounts = sum(data_use, 1);
    weightedSum = sum(bsxfun(@times, data_use, LTTime_use), 1);
    
    valid = totalCounts > 0;
    MPET(valid) = weightedSum(valid) ./ totalCounts(valid) - t0;

end




%%%%%%%%%%%%%%%%%%%%

function [totalBG, effectiveBG, afterpulseBG, correctedData] = internal_calcEffectiveBG(LTHistgram, measuredBG, deadTime, samplingFreq, afterpulseRatio)

% this is to standardize calculation of background with two considerations:
% 1. After pulse ratio
% 2. the background is smaller at high intensity due to deadtime.
% Note, measureBG is nx1 vector. LTHistgram is nxmxch, n is LT times, m is
% sampling time, and ch are number of channels.
% deadtime is in second.

    siz = size(LTHistgram);
    LTHistgram = double(LTHistgram); % its default is uint32.
    
    BGCountEfficiency = 1 - sum(measuredBG) * samplingFreq * deadTime; % E = 1/(1+r0.*dt) = 1-r.*dt;
    trueBG = measuredBG ./ BGCountEfficiency;
    
    rawIntensityPerSecond = sum(LTHistgram, 1) * samplingFreq;
    countEfficiency = 1 - rawIntensityPerSecond.*deadTime; % E = 1/(1+r0.*dt) = 1-r.*dt;
    
    effectiveBG = repmat(trueBG, [1, siz(2), 1]) .* repmat(countEfficiency, [siz(1), 1, 1]);
    effectiveBGInt = sum(trueBG) .* countEfficiency * samplingFreq;
    
    %note: effective BG includes AP. That is why it needs to be subtracted when
    %calculating afterpulseBG. Otherwise, there will be over calculation.
    % F = F[measured] - BG[measured and corrected for intensity] -
    % APRatio*(Intensity[measured] - BGIntensity[measured and corrected for
    % intensity])/#ofpoints.
    afterpulseBG = repmat((rawIntensityPerSecond - effectiveBGInt) * afterpulseRatio / (samplingFreq * siz(1)), [siz(1), 1, 1]);
    
    totalBG = effectiveBG + afterpulseBG;

    correctedData = LTHistgram - totalBG;

end
