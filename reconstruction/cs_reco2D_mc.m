function [image_out,sense_map] = cs_reco2D_mc(app,kspace_in,nc,Wavelet,TVxy,LR,TVt,TVd,ESPIRiT,CLEAR,SOS)

% app = matlab app
% kspace_in = sorted k-space 
% nc = number of RF receiver coils
% averages_in = number of averages per k-line
% Wavelet = wavelet L1-norm regularization factor
% TVt = total variation in CINE dimension
% TVxy = total variation in xy-dimension regularization
% TVd = total variation in dynamic dimension
% ESPIRiT = reconstruction of multi-coil data with ESPIRiT true or false

% kspace_in = {coil}[CINE, x, y, slice, dynamic]
%                     1    2  3    4       5  
dimc = size(kspace_in{1},1);
dimx = 2^nextpow2(size(kspace_in{1},2));
dimy = 2^nextpow2(size(kspace_in{1},3));
dimz = size(kspace_in{1},4);
dimd = size(kspace_in{1},5);

% for convenience make rectangular matrix size
mdimxy = max([dimx,dimy]);
dimx = mdimxy;
dimy = mdimxy;

% resize k-space to next power of 2
for i = 1:nc 
    kspace_in{i} = bart(['resize -c 1 ',num2str(dimx),' 2 ',num2str(dimy)],kspace_in{i});
end

% put all data in a normal matrix
for l = 1:dimd
    for k = 1:dimz
        for j = 1:dimc
            for i = 1:nc
                kspace(j,:,:,k,l,i) = kspace_in{i}(j,:,:,k,l);
            end
        end
    end
end

% Bart dimensions  Bart   Matlab
% 	READ_DIM,       0       1   z  
% 	PHS1_DIM,       1       2   y  
% 	PHS2_DIM,       2       3   x  
% 	COIL_DIM,       3       4   coils
% 	MAPS_DIM,       4       5   sense maps
% 	TE_DIM,         5       6
% 	COEFF_DIM,      6       7
% 	COEFF2_DIM,     7       8
% 	ITER_DIM,       8       9
% 	CSHIFT_DIM,     9       10
% 	TIME_DIM,       10      11  cardiac / respiratory CINE frames
% 	TIME2_DIM,      11      12  dynamics
% 	LEVEL_DIM,      12      13
% 	SLICE_DIM,      13      14  slices
% 	AVG_DIM,        14      15

kspace_pics = permute(kspace,[7,3,2,6,8,9,10,11,12,13,1,5,14,4]);


% wavelet in spatial dimensions 2^1+2^2=6
% total variation in spatial dimensions 2^1+2^2=6
% total variation in cine dimension 2^10 = 1024
% total variation in dynamic dimension 2^11 = 2048

if ESPIRiT && nc>1
    
    % ESPIRiT reconstruction
    TextMessage(app,'ESPIRiT reconstruction ...');
    
    % Calculate coil sensitivity maps with ecalib bart function
    kspace_pics_sum = sum(kspace_pics,[11,12]);
    sensitivities = bart('ecalib -S -I -a', kspace_pics_sum);      % ecalib with softsense
    
    picscommand = 'pics -S ';
    if Wavelet>0
       picscommand = [picscommand,' -RW:6:0:',num2str(Wavelet)];
    end
    if TVxy>0
       picscommand = [picscommand,' -RT:6:0:',num2str(TVxy)];
    end
    if LR>0
       picscommand = [picscommand,' -RL:6:6:',num2str(LR)];
    end
    if TVt>0
       picscommand = [picscommand,' -RT:1024:0:',num2str(TVt)];
    end
    if TVd>0
       picscommand = [picscommand,' -RT:2048:0:',num2str(TVd)];
    end
    image_reg = bart(picscommand,kspace_pics,sensitivities);
    
    % Sum of squares reconstruction
    if SOS==1
        image_reg = bart('rss 16', image_reg);
    end
    image_reg = abs(image_reg);
    
    % CLEAR intensity correction
    if CLEAR
        TextMessage(app,'CLEAR correction ...');
        clear_map = sqrt(sum(abs(sensitivities).^2,[4,5]));         % sum of squares sensitivity maps
        data_dims = size(image_reg);                                % size of the bart reconstructed images
        clear_map = repmat(clear_map,[1 1 1 data_dims(4:end)]);     % adjust size of sensitivity maps to size of images
        clear_map(clear_map<0.2) = 0;                               % threshold to avoid division by very low values
        image_reg = image_reg./clear_map;                           % clear corrected image
        image_reg(isnan(image_reg)) = 0;                            % correct for division by zero
        image_reg(isinf(image_reg)) = 0;                            % correct for division by zero
    end
    
else
    
    % Reconstruction without sensitivity correction
    sensitivities = ones(1,dimy,dimx,nc,1,1,1,1,1,1,1,1,1,dimz);
    
    % regular reconstruction
    picscommand = 'pics -S ';
    if Wavelet>0
       picscommand = [picscommand,' -RW:6:0:',num2str(Wavelet)];
    end
    if TVxy>0
       picscommand = [picscommand,' -RT:6:0:',num2str(TVxy)];
    end
    if LR>0
       picscommand = [picscommand,' -RL:6:6:',num2str(LR)];
    end
    if TVt>0
       picscommand = [picscommand,' -RT:1024:0:',num2str(TVt)];
    end
    if TVd>0
       picscommand = [picscommand,' -RT:2048:0:',num2str(TVd)];
    end
    image_reg = bart(picscommand,kspace_pics,sensitivities);
    
    % Take absolute values
    image_reg = abs(image_reg);
    
end

% rearrange to correct orientation: frames, x, y, slices, dynamics, coils
if SOS==1
    
    image_reg = reshape(image_reg,[dimy,dimx,dimc,dimd,dimz]);
    image_out = flip(permute(image_reg,[3,2,1,5,4]),3);
    
else
    
    image_reg = reshape(image_reg,[dimy,dimx,dimc,nc,dimd,dimz]);
    image_out = flip(permute(image_reg,[3,2,1,6,5,4]),3);
    
end

% sense map orientations: x, y, slice, map1, map2
sense_map = flip(permute(abs(sensitivities),[3,2,14,4,5,1,6,7,8,9,10,11,12,13]),2);

% normalize sense map to reasonable value range
sense_map = sense_map*4095/max(sense_map(:));

end