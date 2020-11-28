function traject = zigzag_trajectory(nr_pe)

% zigzag k-space trajectory

for i=1:nr_pe/2
    traj(i)=2*(i-1)+1;
    traj(i+nr_pe/2)=nr_pe-2*i+2;
end

traject=traj;

end