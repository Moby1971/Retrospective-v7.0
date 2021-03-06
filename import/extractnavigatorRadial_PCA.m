function [amplitude,phaseoffset_y] = extractnavigatorRadial_PCA(inputdata,navigatorlength,nr_nav_points,posneg)

% Extracts the navigator data from the raw k-space data
% Outputs a 1D array of doubles


% Size of the input data, 4th dimension is the readout direction which contains the navigator
[nr_rep,dimz,dimy,dimx] = size(inputdata);  


% Determine the phase offset of the individual spokes based on the center navigator point
cnt = 1;
for i = 1:nr_rep
    for j = 1:dimz
        for k = 1:dimy
            phaseoffset(cnt,:) = [k,phase(inputdata(i,j,k,navigatorlength))];
            cnt = cnt + 1;
        end
    end
end
%phaseoffset = sort(phaseoffset,1);
%phaseoffset = [unique(phaseoffset(:,1)),accumarray(phaseoffset(:,1),phaseoffset(:,2),[],@mean)];

phaseoffset_x = phaseoffset(:,1);
phaseoffset_y = phaseoffset(:,2);
%phaseoffset_y = unwrap(phaseoffset_y);
figure(1)
scatter(phaseoffset_x,phaseoffset_y)


% Extract the navigator and put it in a long array
navdata_amplitude = reshape(permute(inputdata,[3,2,1,4]),nr_rep*dimy*dimz,dimx);


if nr_nav_points > 1

    range = round(nr_nav_points/2);
    
    % Take the principal component of the data
    data = navdata_amplitude(:,navigatorlength-range:navigatorlength+range);
    [coeff,~,~] = pca(data);
    dataPCA = data*coeff;
    
    % Take the principal component of the data
    amplitude = abs(dataPCA(:,1))';
        
else
    
    % single nav point
    amplitude = abs(navdata_amplitude(:,navigatorlength))';

end



% Detrend
amplitude(1,:) = detrend(amplitude(1,:));


% Make a guess whether the respiration peaks are positive or negative
nr_elements = length(amplitude);
first_element = round(0.4*nr_elements);
last_element = round(0.6*nr_elements);
max_amplitude = abs(max(amplitude(1,first_element:last_element)));
min_amplitude = abs(min(amplitude(1,first_element:last_element)));
if min_amplitude > max_amplitude
   amplitude = -amplitude; 
end


% Multiple with +1 or -1 depending on switch
amplitude = amplitude*posneg;


end