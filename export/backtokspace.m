function kspace_out = backtokspace(twoD,movie)

% frames, x, y, slices, dynamics
[nf, ~, ~, dimz, nr] = size(movie);
kspace = zeros(size(movie));

if twoD == 1
    
    parfor i = 1:nf
        
        for j = 1:nr
            
            for k = 1:dimz
                
                kspace(i,:,:,k,j) = fft2c_mri(squeeze(movie(i,:,:,k,j)));
                
            end
            
        end
        
    end
    
    % samples, views, views2, slices, echoes (frames), experiments
    kspace_out = permute(kspace,[2,3,6,4,1,5]);
    
else
    
    %movie = circshift(movie,1,4);
    
    parfor i = 1:nf
        
        for j = 1:nr
            
            kspace(i,:,:,:,j) = fft3c_mri(squeeze(movie(i,:,:,:,j)));
            
        end
        
    end
    
    % samples, views, views2, slices, echoes (frames), experiments
    kspace_out = permute(kspace,[2,3,4,6,1,5]);
    
end


end