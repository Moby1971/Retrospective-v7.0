function [image_out,sense_map] = cs_reco3D_mc(app,kspace_in,nc,Wavelet,TVxyz,LR,TVt,TVd,ESPIRiT)

% app = matlab app
% kspace_in = sorted k-space
% nc = number of RF receiver coils
% averages_in = number of averages per k-line
% Wavelet = wavelet L1-norm regularization factor
% TVt = total variation in time regularization
% TVxyz = total variation in xyz-dimension regularization
% VNSA = variable number of signal averages correction true or false
% ESPIRiT = reconstruction of multi-coil data with ESPIRiT true or false



% kspace_in = {coil}[CINE, x, y, slice, dynamic]
%                     1    2  3    4       5  
dimc = size(kspace_in{1},1);
dimx = 2^nextpow2(size(kspace_in{1},2));
dimy = 2^nextpow2(size(kspace_in{1},3));
dimz = 2^nextpow2(size(kspace_in{1},4));
dimd = size(kspace_in{1},5);

% for convenience make rectangular matrix size
mdimxy = max([dimx,dimy]);
dimx = mdimxy;
dimy = mdimxy;

% resize to next power of 2
for i = 1:nc
    kspace_in{i} = bart(['resize -c 1 ',num2str(dimx),' 2 ',num2str(dimy),' 3 ',num2str(dimz)],kspace_in{i});
end


% bart (kz,ky,kx,channels,multiple sense, ..., ..., cardiac frames at pos 11)
%for i = 1:nc
%    kspace(:,:,:,:,:,i) = kspace_in{i};
%end

for l = 1:dimd
    for k = 1:dimz
        for j = 1:dimc
            for i = 1:nc
                kspace(j,:,:,k,l,i) = kspace_in{i}(j,:,:,k,l);
            end
        end
    end
end



% Bart dimensions
% 	READ_DIM,       1   z  
% 	PHS1_DIM,       2   y  
% 	PHS2_DIM,       3   x  
% 	COIL_DIM,       4   coils
% 	MAPS_DIM,       5   sense maps
% 	TE_DIM,         6
% 	COEFF_DIM,      7
% 	COEFF2_DIM,     8
% 	ITER_DIM,       9
% 	CSHIFT_DIM,     10
% 	TIME_DIM,       11  cardiac / respiratory CINE
% 	TIME2_DIM,      12  dynamics
% 	LEVEL_DIM,      13
% 	SLICE_DIM,      14  slices
% 	AVG_DIM,        15

kspace_pics = permute(kspace,[4,3,2,6,7,8,9,10,11,12,1,5,13,14]);


% wavelet in spatial dimensions 2^0+2^1+2^2=7
% total variation in spatial dimensions 2^0+2^1+2^2=7
% total variation in time 2^10 = 1024
% total variation in dynamic dimension 2^11 = 2048


if ESPIRiT && nc>1
    
    TextMessage(app,'ESPIRiT reconstruction ...');
    kspace_pics_sum = sum(kspace_pics,[11,12]);
    sensitivities = bart('ecalib -I -S -a', kspace_pics_sum);
    
    app.ProgressGauge.Value = 25;
    drawnow;
    
    picscommand = 'pics -S ';
    if Wavelet>0
       picscommand = [picscommand,' -RW:7:0:',num2str(Wavelet)];
    end
    if TVxyz>0
       picscommand = [picscommand,' -RT:7:0:',num2str(TVxyz)];
    end
    if LR>0
       picscommand = [picscommand,' -RL:7:7:',num2str(LR)];
    end
    if TVt>0
       picscommand = [picscommand,' -RT:1024:0:',num2str(TVt)];
    end
    if TVd>0
       picscommand = [picscommand,' -RT:2048:0:',num2str(TVd)];
    end
    image_reg = bart(picscommand,kspace_pics,sensitivities);
    
    app.ProgressGauge.Value = 95;
    drawnow;
    image_reg = bart('rss 16', image_reg);
    image_reg = abs(image_reg);
    
else
    
    % Reconstruction without sensitivity correction
    sensitivities = ones(dimz,dimy,dimx,nc,1,1,1,1,1,1,1,1,1,1);
  
    picscommand = 'pics -S ';
    if Wavelet>0
       picscommand = [picscommand,' -RW:7:0:',num2str(Wavelet)];
    end
    if TVxyz>0
       picscommand = [picscommand,' -RT:7:0:',num2str(TVxyz)];
    end
    if LR>0
       picscommand = [picscommand,' -RL:7:7:',num2str(LR)];
    end
    if TVt>0
       picscommand = [picscommand,' -RT:1024:0:',num2str(TVt)];
    end
    if TVd>0
       picscommand = [picscommand,' -RT:2048:0:',num2str(TVd)];
    end
    
    image_reg = abs(bart(picscommand,kspace_pics,sensitivities));
    
end


% rearrange to correct orientation: frames, x, y, z, dynamics
image_reg = reshape(image_reg,[dimz,dimy,dimx,dimc,dimd]);
image_out = flip(permute(image_reg,[4,3,2,1,5]),3);


% sense map orientations: x, y, z, map1, map2
sense_map = flip(permute(abs(sensitivities),[3,2,1,4,5,6,7,8,9,10,11,12,13,14]),2);


% normalize sense map to reasonable value range
sense_map = sense_map*4095/max(sense_map(:));



app.ProgressGauge.Value = 100;
drawnow;

end