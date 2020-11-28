function export_dicom(directory,movie,parameters,heart_rate,acq_dur,tag)



% export the dicom images

nr_frames = size(movie,1);
dimx = size(movie,2);
dimy = size(movie,3);
dimz = size(movie,4);
nr_dynamics = size(movie,5);

dcmid = dicomuid;   % unique identifier
dcmid = dcmid(1:50);


% create folder if not exist, and clear

folder_name = [directory,[filesep,'CINE_DICOM_',num2str(nr_frames),'_',num2str(dimz),'_',num2str(nr_dynamics),'_',tag]];
if (~exist(folder_name, 'dir')); mkdir(folder_name); end
delete([folder_name,filesep,'*']);


cnt = 1;

for dyn = 1:nr_dynamics
    
    dynamic = ['00',num2str(dyn)];
    dynamic = dynamic(end-2:end);
    
    for slice = 1:dimz
        
        for frame = 1:nr_frames
            
            dcm_header = generate_dicomheader(parameters,slice,nr_dynamics,frame,nr_frames,dimz,dimy,dimx,heart_rate,acq_dur,dcmid);
            fn = ['0000',num2str(cnt)];
            fn = fn(size(fn,2)-4:size(fn,2));
            fname = [folder_name,filesep,'frame_',num2str(frame),'_slice_',num2str(slice),'_dynamic_',num2str(dyn),'_',fn,'.dcm'];
            dicomwrite(squeeze(cast(round(movie(frame,:,:,slice,dyn)),'uint16')), fname, dcm_header);
            
            cnt = cnt + 1;
            
        end
        
    end
    
end


end