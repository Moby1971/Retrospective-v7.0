function [peaks,locations] = trigpoints_cardiac(inputdata,tr,rate,interpolationfactor)

% extracts the ECG/resp trigger points from the navigators
% tr = repetition time TR in ms
% rate = expected heart or respiration rate in bpm
% factor = interpolation factor to make sure that trigger points are not
% exact mulitples of the repetition time

% minimal distance 50% of expected heart rate [in points]   
dist = 0.50*(60/rate)/(tr/1000);      

% cubic spline interpolation of the data
nrl = length(inputdata);
navi = interp1(1:nrl,inputdata(1:nrl),1:1/interpolationfactor:nrl,'spline');

% find the peaks and locations
[peaks,locs]=findpeaks(navi,'MinPeakDistance',dist*interpolationfactor);
locs = locs + interpolationfactor/2;

% recalculate orginal fractional time point peak positions
locations = locs/interpolationfactor;

% locations are in units of samples (actual time is thus locations*TR)


end
