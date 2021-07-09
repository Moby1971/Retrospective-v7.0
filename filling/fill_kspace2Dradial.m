function [kspace,averages,mean_navigator]=fill_kspace2Dradial(raw,includewindow,nr_of_card_frames,nr_of_resp_frames,nr_dynamics,dimz,dimy,nr_ksteps,dimx,card_bin_ass,resp_bin_ass,trajectory,phaseoffset)

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

% dynamics assignment
totalk = nr_reps * nr_ksteps * dimz;
dyn_bin_ass = round(linspace(0.5, nr_dynamics+0.49, totalk));       % list of increasing integer number 1 .. nr_dynamics evenly spaced over the entire acquistion time


disp(trajectory)

cnt = 0;

for slice=1:dimz                    % loop over slices
    
    for r=1:nr_reps                   % loop through all repetitions
        
        for y=1:nr_ksteps            % loop through all the phase-encoding steps
            
            cnt = cnt + 1;
            
            if (card_bin_ass(cnt) > 0) && (includewindow(cnt) == 1)      % if assigment = 0, this acquisition is discarded
                
                kline = trajectory(mod(cnt - 1,nr_ksteps) + 1);    % the phase-encoding step using the trajectory info
                sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) = sorted_kspace(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) + unsorted_kspace(1,r,slice,y,:,1)*exp(-1i*phaseoffset(cnt));   % add the data to the correct k-position
                sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) = sorted_averages(resp_bin_ass(cnt),card_bin_ass(cnt),slice,kline,:,dyn_bin_ass(cnt)) + 1;        % increase the number of averages with 1
                                                
            end
            
        end
        
    end
    
end


% find center of k-space
kspacesum = squeeze(sum(sorted_kspace,[1 2 3 6]));                 % sum over all slices frames and dynamics
[row, col] = find(ismember(kspacesum, max(kspacesum(:))));         % coordinate of k-space maximum = center of k-space


sorted_kspace = sorted_kspace./sorted_averages;   % normalize by number of averages
sorted_kspace(isnan(sorted_kspace)) = complex(0); % correct for NaN or Inf because of division by zero in case of missing k-lines
sorted_kspace(isinf(sorted_kspace)) = complex(0);


disp(size(sorted_kspace));

imshow(squeeze(real(sorted_kspace(1,1,1,:,:))),[]);

%plot(squeeze(abs(kspace(1,1,1,:,64))));

% Apply a circular Tukey filter
filterwidth = 0.2;
flt = circtukey2D(dimy,dimx,row,col,filterwidth);
tukeyfilter(1,1,1,:,:,1) = flt;

kspace = sorted_kspace.*tukeyfilter;
averages = sorted_averages;

end