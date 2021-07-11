function [kspace,averages] = fill_kspace2D(app,raw,includewindow,nr_of_card_frames,nr_of_resp_frames,nr_dynamics,dimz,dimy,nr_ksteps,dimx,card_bin_ass,resp_bin_ass,trajectory,share)

% This function creates 3 arrays
% (1) the kspace data sorted into the correct cardiac frames and phase-encoding positions
% (2) an array with the same size that keeps track of the number of averages per k-space point for normalization/statistics purposes
% (3) an average cardiac navigator signal

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

sorted_kspace = complex(zeros(nr_of_resp_frames,nr_of_card_frames,dimz,dimy,dimx,nr_dynamics));    % fill temp k-space with zeros
sorted_averages = zeros(nr_of_resp_frames,nr_of_card_frames,dimz,dimy,dimx,nr_dynamics);  % fill temp nr averages array with zeros
nr_reps = size(raw,1);                           % number of k-space repetitions
unsorted_kspace = reshape(raw,[1,size(raw),1]);


% Dynamics assignment
totalk = nr_reps * nr_ksteps * dimz;
dyn_bin_ass = round(linspace(0.5, nr_dynamics+0.49, totalk));       % list of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time



% Sorting
cnt = 0;

for slice=1:dimz                    % loop over slices
    
    for i=1:nr_reps                 % loop through all repetitions
        
        for j=1:nr_ksteps           % loop through all the phase-encoding steps
            
            cnt = cnt + 1;
            
            if (card_bin_ass(cnt) > 0) && (includewindow(cnt) == 1)      % if assigment = 0, this acquisition is discarded
                
                kline = trajectory(mod(cnt - 1,nr_ksteps) + 1);    % the phase-encoding step using the trajectory info
                sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) = sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) + unsorted_kspace(1,i,slice,j,:,1);   % add the data to the correct k-position
                sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) = sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) + 1;        % increase the number of averages with 1
                
            end
            
        end
        
    end
    
end



% Temp k-space
new_kspace = sorted_kspace;
new_averages = sorted_averages;



% Find center of k-space
kspacesum = squeeze(sum(sorted_kspace,[1 2 3 6]));                 % sum over all slices frames and dynamics
[row, col] = find(ismember(kspacesum, max(kspacesum(:))));         % coordinate of k-space maximum = center of k-space



% Weighted view sharing
if (share > 0) && (nr_of_card_frames > 1 || nr_of_resp_frames > 1) 
    
    app.TextMessage('View sharing ...');
    
    % respiratory or cardiac frames
    nrframes = nr_of_card_frames;
    if nr_of_resp_frames > 1
       nrframes = nr_of_resp_frames;
       new_kspace = permute(new_kspace,[2,1,3,4,5,6]);
       new_averages = permute(new_averages,[2,1,3,4,5,6]);
       sorted_kspace = permute(sorted_kspace,[2,1,3,4,5,6]);
       sorted_averages = permute(sorted_averages,[2,1,3,4,5,6]);
    end
    
    % determine share range
    maxshare = round(max([nr_of_card_frames nr_of_resp_frames])/2);     % maximum number of shares
    share(share > maxshare) = maxshare;
        
    % define ellipsoid regions
    Ry = round(dimy/share/2);
    Rx = round(dimx/share/2);
    [Y,X] = ndgrid(1:dimy,1:dimx);
    for i = 1:share
        L(i,:,:) = sqrt( ((row-Y)/(Ry*i)).^2 + ((col-X)/(Rx*i)).^2 ) <= 1;
    end
    C(1,:,:) = L(1,:,:);
    if share > 1
        for i = 2:share
            C(i,:,:) = L(i,:,:) - L(i-1,:,:);
        end
    end
    
    % weights
    for i = 1:share
        for j = 1:share
            weights(i,j) = gauss(i+j-1,share,0);
        end
    end
    weights = 0.5*weights/max(weights(:));
    
    disp(weights)
    
    % apply sharing to k-space
    for frame = 1:nrframes
        
        for i = -share:share
            
            sharedframe = frame + i;
            sharedframe(sharedframe < 1) = nrframes - sharedframe - 1;
            sharedframe(sharedframe > nrframes) = sharedframe - nrframes;
            
            if i~=0
                
                for j = 1:share
                    
                    ROI = reshape(squeeze(C(j,:,:)),[1 1 1 dimy dimx 1])*weights(j,abs(i));
                    new_kspace(:,frame,:,:,:,:)   = new_kspace(:,frame,:,:,:,:)   + sorted_kspace(:,sharedframe,:,:,:,:)   .* ROI;
                    new_averages(:,frame,:,:,:,:) = new_averages(:,frame,:,:,:,:) + sorted_averages(:,sharedframe,:,:,:,:) .* ROI;
                    
                end
                
            end
            
        end
        
    end
    
    % respiratory or cardiac frames
    if nr_of_resp_frames > 1
       new_kspace = permute(new_kspace,[2,1,3,4,5,6]);
       new_averages = permute(new_averages,[2,1,3,4,5,6]);
    end

end



% Normalize by number of averages
new_kspace = new_kspace./new_averages;   
new_kspace(isnan(new_kspace)) = complex(0); % correct for NaN or Inf because of division by zero in case of missing k-lines
new_kspace(isinf(new_kspace)) = complex(0);



% Apply a circular Tukey filter
filterwidth = 0.1;
flt = circtukey2D(dimy,dimx,row,col,filterwidth);
tukeyfilter(1,1,1,:,:,1) = flt;



% Report back
kspace = new_kspace.*tukeyfilter;
averages = new_averages;


end