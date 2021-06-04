function export_dicom_dcm(app,dcmdir,directory,movie,parameters,heart_rate,acq_dur,tag,recotype)



% size of the data
nr_frames = size(movie,1);
dimx = size(movie,2);
dimy = size(movie,3);
dimz = size(movie,4);
nr_dynamics = size(movie,5);


% create folder if not exist, and delete folder content
folder_name = [directory,[filesep,'CINE_DICOM_',num2str(nr_frames),'_',num2str(dimz),'_',num2str(nr_dynamics),'_',tag]];
if (~exist(folder_name, 'dir')); mkdir(folder_name); end
delete([folder_name,filesep,'*']);


% reading in the DICOM header information
listing = dir(fullfile(dcmdir, '*.dcm'));
dcmfilename = [listing(1).folder,filesep,listing(1).name];
base_header = dicominfo(dcmfilename);
TextMessage(app,strcat('Reading DICOM info from',{' '},dcmfilename));

% export the dicom images

cnt = 1;


for dyn = 1:nr_dynamics
    
    dynamic = ['00',num2str(dyn)];
    dynamic = dynamic(end-2:end);
    
    for slice = 1:dimz
        
        for frame=1:nr_frames
            
            % read header for specific slice location
            loc = dimz+1-slice;
            dcmfilename = [listing(loc).folder,filesep,listing(loc).name];
            base_header = dicominfo(dcmfilename);
            
            % generate a dicom header
            dcm_header = generate_dicomheader_dcm(base_header,parameters,nr_dynamics,frame,slice,dyn,nr_frames,dimz,dimy,dimx,heart_rate,acq_dur,recotype);
            % filename
            fn = ['0000',num2str(cnt)];
            fn = fn(size(fn,2)-4:size(fn,2));
            fname = [folder_name,filesep,'frame_',num2str(frame),'_slice_',num2str(slice),'_dynamic_',num2str(dyn),'_',fn,'.dcm'];
            
            % write dicom file
            dicomwrite(squeeze(cast(round(movie(frame,:,:,slice,dyn)),'uint16')), fname, dcm_header);
            
            cnt = cnt + 1;
            
        end
        
    end
    
end



end