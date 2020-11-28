function [card_bin_nrs,resp_bin_nrs] = assign_bin_frames_abs(bin_times_card,bin_times_resp,resp_window,nr_klines,nr_resp_frames)

% assigns all the measured k-lines to a specific cardiac phase and respiratory phase bin
% 
% ABSOLUTE BINNING

card_assignments = zeros(nr_klines,1);                           % zero = no assignment (breathing, begin/end of data)
resp_assignments = zeros(nr_klines,1);                           % zero = no assignment (breathing, begin/end of data)
startpoint = round(bin_times_resp(1,2)+1);                       % first measurement to be considered
endpoint = round(bin_times_resp(1,size(bin_times_resp,2))-1);    % last measurement to be considered
loc_maxj_card=size(bin_times_card,2);                            % last bin border
loc_maxj_resp=size(bin_times_resp,2);

% Assignmnent of all k-line acquisitions to a specific bin
%
% INPUT : bin_times_card = array of all time-stamps (units of samples) of the whole acquisition
% 
%                       >    <--- decreases ------ bin_times_card(j)    
%                      ||            |
%                      ||  j-th bin  |
%                      ||            |
%                      ||            |
%       - increases -----> i-th time point / sample
%
% OUTPUT: For each time point / k-line sample assignment to cardiac bin (saw-tooth pattern)
%
%
%
%         /|      /|    /|        /|           
%        / |     / |   / |      /  |
%       /  |   /   |  /  |    /    |  /     = ABSOLUTE BINNING, time between frames is equal
%      /   |  /    | /   |  /      | /                          despite differences in heart beat duration 
%     /    |/      |/    |/        |/                           -> later frames receive less data
%     12345 1234567 12345 123456789  <- frames/bins
%


parfor i=startpoint:endpoint                                     % start search to which heartbeat the measurement belongs
    j=loc_maxj_card;
   
    while i < bin_times_card(1,j)
        j=j-1;
    end
    
    card_assignments(i) = bin_times_card(2,j);                   % assign to bin frame number 
   
    if nr_resp_frames==1
        if resp_window(i)==1                                     % if measurement is during respiration and only 1 resp state, put back to 0 to discard this k-line
            card_assignments(i) = 0; 
        end
    end
    
end

% Respiratory binning is still done using PHASE BINNING

parfor i=startpoint:endpoint                                     % start search to which respiration the measurement belongs
    j=loc_maxj_resp;
    while i<bin_times_resp(1,j)
        j=j-1;
    end
    resp_assignments(i) = mod(j-1,nr_resp_frames)+1;             % assign to bin frame number = j modulus nr_frames
end


card_bin_nrs = card_assignments;
resp_bin_nrs = resp_assignments;

end