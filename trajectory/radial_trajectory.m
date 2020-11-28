function traject = radial_trajectory(nr_pe)

% Radial k-space trajectory
% The actual trajectory is later incorporated in the trajectory that is used by the Bart toolbox
% Therefore the trajectory is simply linear here

for i=1:nr_pe
    traject(i) = i;
end

end