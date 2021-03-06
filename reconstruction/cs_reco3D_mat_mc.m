function image_out = cs_reco3D_mat_mc(app,kspace_in,nc,averages,lambda_TV)


% kdata_in = {coil}[CINE, y, x, z, dynamic]
%                     1   2  3  4     5
nr_cine = size(kspace_in{1},1);
nr_dynamics = size(kspace_in{1},5);


% in case of 1 frame, duplicate that frame to facilitate reconstruction
if nr_cine == 1
    for i = 1:nc
        kspace_in{i}(2,:,:,:,:) = kspace_in{i}(1,:,:,:,:);
    end
end


% kspace data z,y,x,cine,dynamic
for i = 1:nc
    kspace_in{i} = permute(kspace_in{i},[4,2,3,1,5]);
end


% kspace data z,y,x,cine,dynamic,coils
for i = 1:nc
    kspace(:,:,:,:,:,i) = kspace_in{i};
end


% averages data z,y,x,cine,dynamics
averages = permute(averages,[4,2,3,1,5]);


% reset progress counter
param.iteration = 0;


for dynamic = 1:nr_dynamics
    
    % kspace of slice and dynamic
    kdata = squeeze(kspace(:,:,:,:,dynamic,:));
 
    % pad to next power of 2
    dimz = 2^nextpow2(size(kdata,1));
    dimy = 2^nextpow2(size(kdata,2));
    dimx = 2^nextpow2(size(kdata,3));
    
    % for convenience make rectangular matrix size
    mdimxy = max([dimx,dimy]);
    dimx = mdimxy;
    dimy = mdimxy;

    padsizez = round((dimz - size(kdata,1))/2);
    padsizey = round((dimy - size(kdata,2))/2);
    padsizex = round((dimx - size(kdata,3))/2);
    kdata = padarray(kdata,[padsizez,padsizey,padsizex,0],'both');
    kdata = kdata(1:dimz,1:dimy,1:dimx,:,:);
    
    % size of the data, z,y,x,frames,coils
    [nz,ny,nx,~,nc] = size(kdata);
    
    % normalize the data in the range of approx 0 - 1 for better numerical stability
    kdata = kdata/max(abs(kdata(:)));
    
    % kspace mask, 0 = nodata, 1 = data, zero-pad to same size as k-space
    mask = squeeze(averages(:,:,:,:,dynamic));
    mask = padarray(mask,[padsizez,padsizey,padsizex,0],'both');
    mask = mask(1:dimz,1:dimy,1:dimx,:);
    mask = mask./mask;
    mask(isnan(mask)) = 1;
    mask = logical(mask);
    
    % coil sensitivity map
    b1 = ones(nz,ny,nx,nc);
    
    % data
    param.y = kdata;
    
    % reconstruction design matrix
    param.E = Emat_zyxt(mask,b1);
    
    % Total variation (TV) constraint in the temporal domain
    % for 'consistency' with Bart reconstruction, TV seems to be scale empirically by a factor of 8
    % TV only in the time domain
    param.TV = TVOP3D;
    param.TVWeight = lambda_TV/8;
    
    % number of iterations, 2 x 10 iterations
    param.nite = 10;
    param.nouter = 2;
    param.totaliterations = nr_dynamics * param.nouter * param.nite;
    
    % linear reconstruction
    kdata1 = randn(size(kdata))/2000 + kdata;  % add a little bit of randomness, such that linear reco is not exactly right
    recon_dft = param.E'*kdata1;
    
    % iterative reconstruction
    recon_cs=recon_dft;
    for n=1:param.nouter
        [recon_cs,param.iteration] = CSL1NlCg(app,recon_cs,param);
    end
    
    % rearrange to correct orientation: frames, x, y, z
    image_tmp = flip(permute(abs(recon_cs),[4, 2, 3, 1]),2);

    % output reconstructed image
    image_out(:,:,:,:,dynamic) = image_tmp;
    
end


% correct back to 1 frame reconstruction
if nr_cine == 1
   image_out = image_out(1,:,:,:,:);
end


% update gauge
app.ProgressGauge.Value = 100;
drawnow;


end