function [peaks,locations] = trigpoints_resp(inputdata,~,~,interpolationfactor)

% extracts the ECG/resp trigger points from the navigators
% tr = repetition time TR in ms
% rate = expected heart or respiration rate in bpm
% factor = interpolation factor to make sure that trigger points are not
% exact mulitples of the repetition time

% minimal distance 55% of expected rate [in points] to prevent 2nd harmonic detection   
%dist = 0.55*(60/rate)/(tr/1000);      

% cubic spline interpolation of the data
nrl = size(inputdata,1);
navi = interp1(1:nrl,inputdata(1:nrl),1:1/interpolationfactor:nrl,'spline');

% find the peaks and locations
[peaks,locs]=findpeaks(navi,'MinPeakProminence',0.1);
locs = locs + interpolationfactor/2;

% recalculate orginal fractional time point peak positions
locations = locs/interpolationfactor;

end
