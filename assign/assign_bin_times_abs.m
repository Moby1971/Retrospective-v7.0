function [card_bins,resp_bins]=assign_bin_times_abs(card_locations,nr_of_card_frames,resp_locations,nr_of_resp_frames,heartrate,TR)

nr_card = size(card_locations,2); % number of cardiac time points
nr_resp = size(resp_locations,2); % number of respiratory time points


%  cardiac binning = ABSOLUTE BINNING
%
%  INPUT: card_locations = array with the cardiac trigger points in units of samples
%         nr_of_card_frames = number of desired cardiac bins
%
%  card_locations(i)                                card_locations(i+1)
%       ||  j=1     j=2     j=3     j=4     j=5     j=nr_of_card_frames
%       ||       |       |       |       |       |       ||
%       ||       |       |       |       |       |       ||
%       ||       |       |       |       |       |       ||
%       cnt=1   cnt=2   cnt=3   cnt=4   cnt=5   cnt=6   cnt=7   cnt=...
%      cbins(1) cbins(2) .....
%       
%  RESULT: cbins = array with time-stamp of all the cardiac bins for all heartbeats in the measurement in units of samples
%

% frame duration in sample tijd
framedur = (60/heartrate)*(1000/TR)/nr_of_card_frames;

cbins = zeros(2,(nr_card-1)*nr_of_card_frames);

cnt = 1;
for i=1:nr_card-1

    nrframes = floor((card_locations(i+1)-card_locations(i))/framedur); % number of frames for this particular heartbeat
    if nrframes>nr_of_card_frames 
        nrframes = nr_of_card_frames; 
    end
    
    for j=1:nrframes
        cbins(1,cnt) = card_locations(i)+(j-1)*framedur; % divide heart-beat in nrframes equal bins of duration framedur
        cbins(2,cnt) = j;                                % immediately assign frame number
        cnt = cnt + 1;
    end
    
end

cbins = cbins(:,1:cnt-1);



% respiratory binning = PHASE BINNING

rbins = zeros(1,(nr_resp-1)*nr_of_resp_frames);

cnt = 1;
for i=1:nr_resp-1
  for j=1:nr_of_resp_frames
       rbins(cnt)=resp_locations(i)+(j-1)*(resp_locations(i+1)-resp_locations(i))/nr_of_resp_frames;
       cnt = cnt + 1;
  end
end


card_bins = cbins;
resp_bins = rbins;

end