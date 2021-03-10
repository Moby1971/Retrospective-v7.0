function [detectedHR,detectedRR,powerspectrum,f] = determine_power_spectrum_PCA(nav_input,nc,TR,dimy,includewindow,excludewindow)

if nc > 1
    
    for i = 1:nc
        data(:,i) = nav_input{i};
    end
    
    % Take the principal component of the data
    [coeff,~,~] = pca(data);
    dataPCA = data*coeff;
    nav_amplitude = dataPCA(:,1);
  
else
    
    nav_amplitude = nav_input{1}';
        
end

% include only those navigators in when includewindow == 1 and excludewindow == 1
nav_amplitude = nav_amplitude.*includewindow.*excludewindow;

% determine the frequency power spectrum
y = fft(nav_amplitude);
fs = 1000/TR;                               % sample frequency in Hz
n = length(nav_amplitude);                  % number of samples
f = (0:n-1)*(fs/n)*60;                      % frequency range in bpm
power = abs(y).^2/n;                        % power of the DFT

% determine frequency and harmonics of k-space trajectory and set those to zero
kfreq = 60/(0.001*dimy*TR);
ifreq = (fs/n)*60;

for i = 1:10
   power(round(i*kfreq/ifreq)+1) = 0;
   power(round(i*kfreq/ifreq)) = 0;
   power(round(i*kfreq/ifreq)-1) = 0;
end

power = movmean(power,3);   % smooth the power spectrum with moving average
powerspectrum = power;

% detect heart rate
minheartbpm = 350;
maxheartbpm = 600;
minidx = round(minheartbpm*n/(fs*60));
maxidx = round(maxheartbpm*n/(fs*60));
[~, idx] = max(power(minidx:maxidx));
detectedHR = round(idx*fs*60/n + minheartbpm);

% detect respiratory rate
minRRbpm = 30;
maxRRbpm = 70;
minidx = round(minRRbpm*n/(fs*60));
maxidx = round(maxRRbpm*n/(fs*60));
[~, idx] = max(power(minidx:maxidx));
detectedRR = round(idx*fs*60/n + minRRbpm);

end