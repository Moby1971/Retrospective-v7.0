function [heart_outputdata,resp_outputdata,pc_outputdata] = filter_nav_PCA(app,inputdata,ncoils,TR,hr,bwhr,rr,bwrr,filtersettings)   

% applies a bandwidth filter on the navigator data
% TR = repetition time TR in ms
% hr = expected heart rate in bpm
% bw = expected bandwidth in heart rate in bpm
% order of the filter

sf=1000/TR;     % sampling frequency in Hz = 1/TR[ms]

resp_harmonics = 2;   % number of higher order harmonics for respiratory frequency, 2 = main + 1 harmonic

order = filtersettings(5); % butterworth filter order

if ncoils > 1
    
    for i = 1:ncoils
        data(:,i) = inputdata{i};
    end
    
    % Take the principal component of the data
    [coeff,~,~] = pca(data);
    dataPCA = data*coeff;
    
    % Take the principal component of the data
    nav_amplitude = dataPCA(:,1);
        
else
    
    nav_amplitude = inputdata{1}';
    
end



% Filter for heart motion
hrf=hr/60;      % expected heartrate in Hz = hr[bpm]/60
bwh=bwhr/60;    % bandwidth heartrate in Hz = [bpm]/60
[b,a] = butter(order,[hrf-0.5*bwh,hrf+0.5*bwh]/(sf/2),'bandpass');   % butterworth bandpass filter
heart_outputdata = filtfilt(b,a,nav_amplitude);   % apply zero-phase shift filtering

% detrend
heart_outputdata = detrend(heart_outputdata);

% normalize envelope
factor = round(5000/app.HeartEditField.Value);  % adapt the envelope setting to match the expected heart rate frequency
[env,~] = envelope(heart_outputdata,factor,'peak');
heart_outputdata = heart_outputdata./abs(env);

            
% Filter for respiration motion
while true
    
    rrf=rr/60;      % expected resprate in Hz = rr[bpm]/60
    bwr=bwrr/60;    % bandwidth resprate in Hz = [bpm]/60
    
    resp_outputdata = zeros(size(nav_amplitude));
    
    if rr<45
        [b, a] = butter(order,(rrf+0.5*bwr)/(sf/2),'low');           % butterworth lowpass filter for low frequencies
        resp_outputdata = filtfilt(b,a,nav_amplitude);
    else
        for i = 1:resp_harmonics
            [b, a] = butter(order,[i*rrf-0.5*bwr,i*rrf+0.5*bwr]/(sf/2),'bandpass');           % butterworth bandpass filter
            resp_outputdata = resp_outputdata + (1/i^2)*filtfilt(b,a,nav_amplitude);          % apply zero-phase shift filtering
        end
    end
    
      
    % In some cases the filter produces NaN when filtering with too low
    % frequency. In those cases the respirate and bandwidth will be
    % increased until there are no more NaN
    
    if sum(isnan(resp_outputdata))==0
        break;
    end
    rr = rr+1;
    app.RespirationEditField.Value = rr;
    bwrr = bwrr+1;
    app.RespirationWidthEditField.Value = bwrr;
    drawnow;
    
end


% detrend and normalize envelope
resp_outputdata = detrend(resp_outputdata);


% normalize envelope
factor = round(5000/app.RespirationEditField.Value);  % adapt the envelope setting to match the expected heart rate frequency
[env,~] = envelope(resp_outputdata,factor,'peak');
resp_outputdata = resp_outputdata./abs(env);

% Return the principal component, or in case 1 coil the orignal navigator data, and detrend
pc_outputdata = detrend(nav_amplitude);


end