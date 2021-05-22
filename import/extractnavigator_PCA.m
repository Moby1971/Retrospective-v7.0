function amplitude = extractnavigator_PCA(inputdata,navigatorlength,nr_nav_points,posneg)

% extracts the navigator data from the raw k-space data
% outputs a 1D array of doubles


% size of the input data, 4th dimension is the readout direction which contains the navigator
[nr_rep,dimz,dimy,dimx] = size(inputdata);   


% extract the navigator and put it in a long array
navdata_amplitude = reshape(permute(inputdata,[3,2,1,4]),nr_rep*dimy*dimz,dimx);



if nr_nav_points > 1

    % Take the principal component of the data
    data = navdata_amplitude(:,navigatorlength-nr_nav_points+1:navigatorlength);
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
last_element = round(0.8*nr_elements);
max_amplitude = abs(max(detrend(amplitude(1,first_element:last_element))));
min_amplitude = abs(min(detrend(amplitude(1,first_element:last_element))));
if min_amplitude > max_amplitude
   amplitude(1,:) = -amplitude(1,:); 
end


% Multiple with +1 or -1 depending on switch
amplitude = amplitude*posneg;


end