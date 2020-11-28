function output = norm_images(input)

% normalize the images to 2^15 range

input = round(32766*input/max(input(:)));

output = input;

end