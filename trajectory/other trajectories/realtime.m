% realtime trajectory test
% Gustav Strijkers, May 2021

clearvars;

% cardiac assignments
z = 1;
for cnt1 = 1:60
    range = 13 + randi(4);
    for cnt2 = 1:range
        hb(z) = cnt1;
        z = z + 1;
    end
end

traj = [];
% kspace trajectory
for i = 1:60
   traj = [traj,zigzag_trajectory(48)];
end


kspace = zeros(60,48,48);

% fill k-space
for i = 1:length(hb)
    
    kspace(hb(i),traj(i),:) = kspace(hb(i),traj(i),:) + 1;

end

kspace2 = zeros(60,48,48);

% sharing
for i = 2:59
    
   kspace2(i,:,:) = kspace(i,:,:) + kspace(i-1,:,:) + kspace(i+1,:,:);
   
end
kspace2(1,:,:) = kspace(1,:,:);
kspace2(60,:,:) = kspace(60,:,:);


figure(1);

for i = 1:60
    subplot(6,10,i),imshow(squeeze(kspace2(i,:,:)),[0 2]);
end


disp(nnz(kspace(:)))
