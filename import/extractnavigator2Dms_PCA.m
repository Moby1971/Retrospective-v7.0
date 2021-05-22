function amplitude = extractnavigator2Dms_PCA(inputdata,navigatorlength,nr_nav_points,posneg)

% extracts the navigator data from the raw k-space data of multi-slice 2D data
% outputs a 1D array of doubles


% size of the input data, 4th dimension is the readout direction which contains the navigator
[nr_rep,dimz,dimy,dimx] = size(inputdata);   


% extract the navigator and put it in a long array
% y-dimension, repetitions, slice, readout
navdata_amplitude = reshape(permute(inputdata,[3,1,2,4]),nr_rep*dimy*dimz,dimx);

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


% Make a guess whether the respiration peaks are positive or negative in the different slices
% This will not be needed with out-of-slice navigator

nr_elements = length(amplitude);
nr_elementsperslice = round(nr_elements/dimz);

for i = 1:dimz

    % First and last navigator point for each slice
    first_element = (i-1)*nr_elementsperslice + 1;
    last_element = i*nr_elementsperslice;
    
    % Only look at part of that data away from the start to prevent transient effects
    first_element1 = (i-1)*nr_elementsperslice + 1 + round(0.4*nr_elementsperslice);
    last_element1 = i*nr_elementsperslice - round(0.1*nr_elementsperslice);
 
    % Min/max of navigator
    max_amplitude = abs(max(detrend(amplitude(1,first_element1:last_element1))));
    min_amplitude = abs(min(detrend(amplitude(1,first_element1:last_element1))));

    if min_amplitude > max_amplitude
        amplitude(1,first_element:last_element) = -amplitude(1,first_element:last_element); 
    end
    
    amplitude(1,first_element:last_element) = detrend(amplitude(1,first_element:last_element));
    
end


% Multiple with +1 or -1 depending on global switch
amplitude = amplitude*posneg;


end