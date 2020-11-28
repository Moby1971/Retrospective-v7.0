function output = calc_resprate(locations,tr)

nrlocs = size(locations,2);

% determine respiratory rate as function of time in bpm

for i=1:nrlocs-1
    rate(i) = 60/((locations(i+1)-locations(i))*(tr/1000));
end

output = rate';

end