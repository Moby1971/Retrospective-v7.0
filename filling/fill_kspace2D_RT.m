function [kspace,averages,nr_dynamics] = fill_kspace2D_RT(app,raw,nr_card_frames,nr_of_resp_frames,dimz,dimy,nr_ksteps,dimx,bin_times_card,trajectory,share)

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
% shared                = weighted view sharing of neighboring data



nr_dynamics = round(length(bin_times_card)/nr_card_frames);                              % number of dynamics equals the number of heartbeats in the acquisition
sorted_kspace = complex(zeros(nr_of_resp_frames,nr_card_frames,dimz,dimy,dimx,nr_dynamics));    % fill temp k-space with zeros
sorted_averages = zeros(nr_of_resp_frames,nr_card_frames,dimz,dimy,dimx,nr_dynamics);           % fill temp averages array with zeros
nr_reps = size(raw,1);                                                                   % number of k-space repetitions
nr_klines = nr_reps * nr_ksteps * dimz;                                                  % total number of k-lines


% unsorted k-lines
raw = permute(raw,[3,1,2,4]);
unsorted_klines = reshape(raw,[1,1,1,nr_klines,dimx]);


% fill k-space
fcnt = 1;
 
for i = 1 : nr_dynamics
    
    for j = 1 : nr_card_frames
     
        for k = 1:nr_klines
            
            if fcnt < length(bin_times_card)
                
                if k > bin_times_card(1,fcnt) && k < bin_times_card(1,fcnt+1)
                    
                    % the phase-encoding step using the trajectory info
                    kline = trajectory(mod(k - 1,nr_ksteps) + 1);
                    
                    % fill k-line
                    sorted_kspace(1,j,1,kline,:,i) = sorted_kspace(1,j,1,kline,:,i) + unsorted_klines(1,1,1,k,:);   % add the data to the correct k-position
                    sorted_averages(1,j,1,kline,:,i) = sorted_averages(1,j,1,kline,:,i) + 1;
                    
                end
                
            end
            
        end
        
        fcnt = fcnt + 1;
                
    end
    
end



% find center of k-space
kspacesum = squeeze(sum(sorted_kspace,[1 2 3 6]));                 % sum over all slices frames and dynamics
[row, col] = find(ismember(kspacesum, max(kspacesum(:))));         % coordinate of k-space maximum = center of k-space



% Temp k-space
new_kspace = sorted_kspace;
new_averages = sorted_averages;



% Weighted view sharing
if (share > 0) && (nr_dynamics > 1) 
    
    app.TextMessage('View sharing ...');
    
    % determine share range
    maxshare = 20;                                          % maximum number of shares
    share(share > maxshare) = maxshare;
    weights = gauss(1:share+1,share,0);
    weights = weights/max(weights);                         % Gaussian weighting
    
    % define ellipsoid regions
    Ry = round(dimy/share/2);
    Rx = round(dimx/share/2);
    [Y,X] = ndgrid(1:dimy,1:dimx);
    for i = 1:share
        L(share - i + 1,:,:) = sqrt( ((row-Y)/(Ry*i)).^2 + ((col-X)/(Rx*i)).^2 ) <= 1;
    end
    
    % apply sharing to k-space
    for frame = 1:nr_dynamics
        
        for i = -share:share
            
            sharedframe = frame + i;
            
            if sharedframe > 0 && sharedframe < nr_dynamics
                
                if i~=0
                    
                    ROI = reshape(squeeze(L(abs(i),:)),[1 1 1 dimy dimx 1])*weights(abs(i));
                    new_kspace(:,:,:,:,:,frame)   = new_kspace(:,:,:,:,:,frame)   + sorted_kspace(:,:,:,:,:,sharedframe)   .* ROI;
                    new_averages(:,:,:,:,:,frame) = new_averages(:,:,:,:,:,frame) + sorted_averages(:,:,:,:,:,sharedframe) .* ROI;
                    
                end
                
            end
            
        end
        
    end

end






% Normalize by number of averages
new_kspace = new_kspace./new_averages;   
new_kspace(isnan(new_kspace)) = complex(0); 
new_kspace(isinf(new_kspace)) = complex(0);



% Apply a circular Tukey filter
filterwidth = 0.1;
flt = circtukey2D(dimy,dimx,row,col,filterwidth);
tukeyfilter(1,1,1,:,:,1) = flt;

kspace = new_kspace.*tukeyfilter;
averages = new_averages;


end