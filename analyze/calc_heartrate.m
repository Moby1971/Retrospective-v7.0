function output = calc_heartrate(locations,tr)

% determine heart rate as function of time in bpm

nrlocs = size(locations,2);

for i=1:nrlocs-1
    rate(i) = 60/((locations(i+1)-locations(i))*(tr/1000));
end

output = rate';

end