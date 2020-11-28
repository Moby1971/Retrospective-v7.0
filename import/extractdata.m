function outputdata = extractdata(inputdata,navigatorlength)

% extracts the imaging k-space data from the raw 2D and 3D k-space data

outputdata = inputdata(:,:,:,navigatorlength+1:end);

end