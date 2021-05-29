function [kspace,averages,nr_dynamics]=fill_kspace2D_RT(raw,nr_card_frames,nr_of_resp_frames,nr_shared_dynamics,dimz,dimy,nr_ksteps,dimx,bin_times_card,trajectory)

% This function creates 3 arrays
% (1) the kspace data sorted into the correct cardiac frames,  and phase-encoding positions
% (2) an array with the same size that keeps track of the number of averages per k-space point for normalization/statistics purposes
% (3) the number of dynamics = number of heart-beats in the datasets

% Required input:
%
% raw                   = unsorted k-space data
% navheart              = navigator signal, used to construct an average heartbeat
% nr_of_card_frames     = number of desired cardiac frames, must be consisted with bin assignments (bin_ass)
% nr_of_resp_frames     = number of desired respiratory frames
% nr_of_dynamics        = number of desired dynamics
% dimz                  = number of slices
% dimy                  = dimensions of the images: dimy (phase encoding)
% nr_ksteps             = number of k-lines in 1 repetition
% dimx                  = dimensions of the images: dimx (readout)
% card_bin_ass          = the cardiac bin assignment array for all measured k-lines
% resp_bin_ass          = the respiratory bin assignment array for all measured k-lines
% traj                  = the k-space trajectory
% includewindow         = data which should be include: 1 = yes, 0 = no



nr_dynamics = round(length(bin_times_card)/nr_card_frames);                              % number of dynamics equals the number of heartbeats in the acquisition
kspace = complex(zeros(nr_of_resp_frames,nr_card_frames,dimz,dimy,dimx,nr_dynamics));    % fill temp k-space with zeros
averages = zeros(nr_of_resp_frames,nr_card_frames,dimz,dimy,dimx,nr_dynamics);           % fill temp averages array with zeros
nr_reps = size(raw,1);                                                                   % number of k-space repetitions
nr_klines = nr_reps * nr_ksteps * dimz;                                                  % total number of k-lines


% unsorted k-lines
raw = permute(raw,[3,1,2,4]);
unsorted_klines = reshape(raw,[1,1,1,nr_klines,dimx]);


% sharing k-lines from neighboring dynamics
shared = [1:nr_shared_dynamics]-round(nr_shared_dynamics/2);


% fill k-space
fcnt = 1;
 
for i = 1 : nr_dynamics
    
    for j = 1 : nr_card_frames
     
        for k = 1:nr_klines
            
            if fcnt < length(bin_times_card)
                
                if k > bin_times_card(1,fcnt) && k < bin_times_card(1,fcnt+1)
                    
                    % the phase-encoding step using the trajectory info
                    kline = trajectory(mod(k - 1,nr_ksteps) + 1);   
                    
                    % share data over dynamics
                    for w = 1 : nr_shared_dynamics
                        dyn = i + shared(w);
                        if dyn > 0 && dyn <= nr_dynamics
                            kspace(1,j,1,kline,:,dyn) = kspace(1,j,1,kline,:,dyn) + unsorted_klines(1,1,1,k,:);   % add the data to the correct k-position
                            averages(1,j,1,kline,:,dyn) = averages(1,j,1,kline,:,dyn) + 1;
                        end
                    end
                    
                end
                
            end
            
        end
        
        fcnt = fcnt + 1;
                
    end
    
end


% Normalize by number of averages
kspace = kspace./averages;   
kspace(isnan(kspace)) = complex(0); 
kspace(isinf(kspace)) = complex(0);


% Apply a circular Tukey filter
filterwidth = 0.2;
flt = circtukey2D(dimy,dimx,filterwidth);
tukeyfilter(1,1,1,:,:,1) = flt;
kspace = kspace.*tukeyfilter;



end