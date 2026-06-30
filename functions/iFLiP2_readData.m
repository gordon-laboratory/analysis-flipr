function Aout = iFLiP2_readData(filename)

% this is to open .iFLiP2 file. 

fid = fopen(filename, 'r');

headerStr = '';
while true
    tline = fgetl(fid);
    assert(ischar(tline), 'end of file reached prematurely');
    if strcmp(tline, 'header_end')
        break;
    end
    headerStr = [headerStr, tline];
end


data = fread(fid, inf, 'uint32=>uint32');
fclose(fid);

try
    evalc(headerStr);
    Aout.header = header;
    nChannels = header.acq.nChannels;

    LTCurveLength = header.acq.LTCurveLength + 1; % extra row is for marks

    dataTimeLength = length(data)/nChannels/LTCurveLength;
    if dataTimeLength~=round(dataTimeLength) % if it is not an integer.
        error('Data size not right!')
    end
    Aout.data = permute(reshape(data, [LTCurveLength, nChannels, dataTimeLength]), [1 3 2]);
    Aout.LTTime = (0:header.acq.LTCurveLength-1)*header.acq.LTResolution;
    Aout.sampleTime = (0.5:dataTimeLength) /  header.state.samplingFreq.Value;
    Aout.marks = Aout.data(end, :, 1); % marks for all channeles are the same.
    Aout.data(end,:,:) = [];

catch
    fprintf('iFLiP2_readData error!');
end

