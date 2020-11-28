function [image_out,sense_map] = cs_reco2D_radial_mc(app,kspace_in,dimx,dimy,ksteps,nc,Wavelet,TVt,TVxy)

% app = matlab app
% kspace_in = sorted k-space 
% nc = number of RF receiver coils
% averages_in = number of averages per k-line
% Wavelet = wavelet L1-norm regularization factor
% TVt = total variation in time regularization
% TVxy = total variation in xy-dimension regularization




% radial k-space trajectory for Bart toolbox

radial = 360;           % 360 degrees
u_spokes = true;        % unique spokes true/false
half_spoke = false;     % half trajectory true/false


for n = 1:dimy
    
    % caluculate angle
    if radial == 180
        a(n) = (n-1)*180/dimy;
    elseif ~u_spokes
        a(n) = (n-1)*360/dimy;
    elseif (n-1) < dimy/2
        a(n) = (n-1)*360/dimy;
    else
        a(n) = (n-1)*360/dimy + 180/dimy;
    end
    
    % calculate x and y values
    for m = 1:dimx
        r(m,n) = (dimx-1)/2 -(m-1);
        x = r(m,n)*-cos(a(n)*(pi/180));
        y = r(m,n)*sin(a(n)*(pi/180));
        traj(1,m,n) = x;
        traj(2,m,n) = y;
        traj(3,m,n) = 0;
    end
    
end


disp(size(traj))

%disp(traj)

% bart (kz,ky,kx,channels,multiple sense, ..., ..., cardiac frames at pos 11)
% sense maps for 2D data:  (kz,ky,kx,....)
for i = 1:nc
   kspace(:,:,:,:,i) = kspace_in{i};    
end

%figure(1)
%imshow(real(squeeze(kspace(1,:,:))),[-2*pi,2*pi]);
%imshow(abs(squeeze(kspace(1,:,:))),[0 10]);

% rearrange for the bart toolbox
kspace_pics = permute(kspace,[4,3,2,5,6,7,8,9,10,11,1]);

% create trajectory
%bartcommand = ['traj -r -y',num2str(ksteps),' -x',num2str(dimx),' -q0:0:0'];
%traj = bart(bartcommand);

disp(size(traj))

%disp(traj)

% sensitivity map
kspace_pics_sum = sum(kspace_pics,11);
lowres_img = bart('nufft -i -l6 -d16:16:1 -t', traj, kspace_pics_sum);
lowres_ksp = bart('fft -u 7', lowres_img);

highres_img = abs(bart('nufft -i -t -l6 -d128:128:1 -t', traj, kspace_pics_sum));
figure(1)
imshow(abs(lowres_img),[0,1.5*max(abs(lowres_img(:)))]);

% zeropad to full size
bartcommand = ['resize -c 0 ',num2str(dimx),' 1 ',num2str(dimx)];
ksp_zerop = bart(bartcommand, lowres_ksp);

% calculate sensitivity map
sensitivities = bart('ecalib -m1', ksp_zerop);

%disp(size(sensitivities))
%sensitivities = ones(128,128);
disp(size(sensitivities))

sensitivities = ones(dimx,dimx,1,1,1,1,1,1,1,1,1,1,1);


picscommand = ['pics -S -u1 -RW:6:0:',num2str(Wavelet),' -RT:6:0:',num2str(TVxy),' -RT:1024:0:',num2str(TVt),' -t'];
image_reg = bart(picscommand,traj,kspace_pics,sensitivities);

image_reg = bart('nufft -i -l1 -d128:128:1 -t', traj, kspace_pics);


disp('*');
disp(size(image_reg));





% rearrange to correct orientation: frames,x, y
image_reg = abs(squeeze(image_reg(:,:,1,1,1,1,1,1,1,1,:)));
image_out = flip(permute(image_reg,[3, 2, 1]),3);
sense_map = flip(permute(squeeze(abs(sensitivities)),[2, 1, 3, 4]),2);



% normalize sense map to reasonable value range
sense_map = sense_map*1024/max(sense_map(:));



end