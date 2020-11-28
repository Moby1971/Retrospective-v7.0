function outputdata = filter_shape(tr,hr,bw,order)   

% determines the filter shape for display purposes
% tr = repetition time TR in ms
% hr = expected heart rate in bpm
% bw = expected bandwidth in heart rate in bpm
% order = order of the filter

sf=1000/tr;    % sampling frequency in Hz = 1/TR[ms]
cf=hr/60;      % expected heartrate in Hz = hr[bpm]/60
bw1=bw/60;      % bandwidth in Hz = [bpm]/60

[b,a] = butter(order,[cf-0.5*bw1,cf+0.5*bw1]/(sf/2),'bandpass');   % butterworth bandpass filter

[h,w] = freqz(b,a,2000,sf);
outputdata = [w*60, abs(h)];

end