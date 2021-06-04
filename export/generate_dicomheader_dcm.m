function dicom_header = generate_dicomheader_dcm(dcmhead,parameters,nr_dynamics,frame,slice,dyn,nr_frames,dimz,dimy,dimx,heart_rate,acq_dur,recotype)

%
% GENERATES DICOM HEADER FOR EXPORT
%
% parameters = parameters from MRD file
% dcmhead = dicom info from scanner generated dicom
%
% frame = current frame number
% nr_frames = total number of frames
% nr_dynamics = number of dynamics
% dimy = y dimension (phase encoding, views)
% dimx = x dimension (readout, samples)
% heart_rate = average heart rat
%
%


if strcmp(recotype,'realtime')
    
    TR = 1000*(60/heart_rate)/nr_frames;    % time between cardiac frames in ms
    TD = 1000*acq_dur/nr_dynamics;                % time between dynamics
    
    pixely = parameters.FOV/dimy;
    pixelx = parameters.FOV/dimx;
    
    fn = ['0000',num2str(frame)];
    fn = fn(size(fn,2)-4:size(fn,2));
    fname = ['cine_',fn,'.dcm'];
    
    dt = datetime(parameters.date,'InputFormat','dd-MMM-yyyy HH:mm:ss');
    year = num2str(dt.Year);
    month = ['0',num2str(dt.Month)]; month = month(end-1:end);
    day = ['0',num2str(dt.Day)]; day = day(end-1:end);
    date = [year,month,day];
    
    hour = ['0',num2str(dt.Hour)]; hour = hour(end-1:end);
    minute = ['0',num2str(dt.Minute)]; minute = minute(end-1:end);
    seconds = ['0',num2str(dt.Second)]; seconds = seconds(end-1:end);
    time = [hour,minute,seconds];
    
    dcmhead.Filename = fname;
    dcmhead.FileModDate = parameters.date;
    dcmhead.FileSize = dimy*dimx*2;
    dcmhead.Width = dimx;
    dcmhead.Height = dimy;
    dcmhead.BitDepth = 15;
    dcmhead.InstitutionName = 'Amsterdam UMC';
    dcmhead.ReferringPhysicianName.FamilyName = 'AMC preclinical MRI';
    dcmhead.StudyDescription = 'CINE';
    dcmhead.InstitutionalDepartmentName = 'Amsterdam UMC preclinical MRI';
    dcmhead.PhysicianOfRecord.FamilyName = 'Amsterdam UMC preclinical MRI';
    dcmhead.PerformingPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
    dcmhead.PhysicianReadingStudy.FamilyName = 'Amsterdam UMC preclinical MRI';
    dcmhead.OperatorName.FamilyName = 'manager';
    dcmhead.ManufacturerModelName = 'MRS7024';
    dcmhead.ReferencedFrameNumber = [];
    %dcmhead.PatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
    %dcmhead.OtherPatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
    dcmhead.RepetitionTime = TR;     % time between frames
    dcmhead.EchoTime = 2;         % approximately correct, unknown because of shortest TE option
    dcmhead.NumberOfAverages = parameters.NO_AVERAGES;
    dcmhead.InversionTime = 0;
    dcmhead.ImagedNucleus = '1H';
    dcmhead.MagneticFieldStrength = 7;
    dcmhead.TriggerTime = (dyn-1)*TD + (frame-1)*TR;    
    dcmhead.AcquisitionMatrix = uint16([dimy 0 0 dimx])';
    dcmhead.AcquisitionDeviceProcessingDescription = '';
    dcmhead.AcquisitionDuration = acq_dur;
    dcmhead.InstanceNumber = (slice-1)*nr_dynamics*nr_frames + (frame-1)*nr_dynamics + dyn;          % instance number
    dcmhead.TemporalPositionIdentifier = (dyn-1)*nr_frames + frame;     
    dcmhead.NumberOfTemporalPositions = nr_dynamics * nr_frames;
    dcmhead.TemporalResolution = TR;
    dcmhead.ImagesInAcquisition = nr_dynamics * nr_frames * dimz;
    dcmhead.TemporalPositionIndex = uint32([]);
    dcmhead.Rows = dimy;
    dcmhead.Columns = dimx;
    dcmhead.PixelSpacing = [pixely pixelx]';
    dcmhead.PixelAspectRatio = [1 1]';
    dcmhead.BitsAllocated = 16;
    dcmhead.BitsStored = 15;
    dcmhead.HighBit = 14;
    dcmhead.PixelRepresentation = 0;
    dcmhead.PixelPaddingValue = 0;
    dcmhead.RescaleIntercept = 0;
    dcmhead.RescaleSlope = 1;
    dcmhead.HeartRate = heart_rate;
    dcmhead.NumberOfSlices = dimz;
    dcmhead.CardiacNumberOfImages = nr_frames;
    
    
    
else
    
    
    
    
    TR = 1000*(60/heart_rate)/nr_frames;    % time between frames in ms
    
    pixely = parameters.FOV/dimy;
    pixelx = parameters.FOV/dimx;
    
    fn = ['0000',num2str(frame)];
    fn = fn(size(fn,2)-4:size(fn,2));
    fname = ['cine_',fn,'.dcm'];
    
    dt = datetime(parameters.date,'InputFormat','dd-MMM-yyyy HH:mm:ss');
    year = num2str(dt.Year);
    month = ['0',num2str(dt.Month)]; month = month(end-1:end);
    day = ['0',num2str(dt.Day)]; day = day(end-1:end);
    date = [year,month,day];
    
    hour = ['0',num2str(dt.Hour)]; hour = hour(end-1:end);
    minute = ['0',num2str(dt.Minute)]; minute = minute(end-1:end);
    seconds = ['0',num2str(dt.Second)]; seconds = seconds(end-1:end);
    time = [hour,minute,seconds];
    
    dcmhead.Filename = fname;
    dcmhead.FileModDate = parameters.date;
    dcmhead.FileSize = dimy*dimx*2;
    dcmhead.Width = dimx;
    dcmhead.Height = dimy;
    dcmhead.BitDepth = 15;
    dcmhead.InstitutionName = 'Amsterdam UMC';
    dcmhead.ReferringPhysicianName.FamilyName = 'AMC preclinical MRI';
    dcmhead.StudyDescription = 'CINE';
    dcmhead.InstitutionalDepartmentName = 'Amsterdam UMC preclinical MRI';
    dcmhead.PhysicianOfRecord.FamilyName = 'Amsterdam UMC preclinical MRI';
    dcmhead.PerformingPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
    dcmhead.PhysicianReadingStudy.FamilyName = 'Amsterdam UMC preclinical MRI';
    dcmhead.OperatorName.FamilyName = 'manager';
    dcmhead.ManufacturerModelName = 'MRS7024';
    dcmhead.ReferencedFrameNumber = [];
    %dcmhead.PatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
    %dcmhead.OtherPatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
    dcmhead.RepetitionTime = TR;     % time between frames
    dcmhead.EchoTime = 2;         % approximately correct, unknown because of shortest TE option
    dcmhead.NumberOfAverages = parameters.NO_AVERAGES;
    dcmhead.InversionTime = 0;
    dcmhead.ImagedNucleus = '1H';
    dcmhead.MagneticFieldStrength = 7;
    dcmhead.TriggerTime = (frame-1)*TR;    % frame time = frame number times calculated TR
    dcmhead.AcquisitionMatrix = uint16([dimy 0 0 dimx])';
    dcmhead.AcquisitionDeviceProcessingDescription = '';
    dcmhead.AcquisitionDuration = acq_dur/nr_dynamics;
    dcmhead.InstanceNumber = (slice-1)*nr_frames + frame;          % instance number
    dcmhead.TemporalPositionIdentifier = frame;     % frame number
    dcmhead.NumberOfTemporalPositions = nr_frames;
    dcmhead.TemporalResolution = TR;
    dcmhead.ImagesInAcquisition = nr_frames*dimz;
    dcmhead.TemporalPositionIndex = uint32([]);
    dcmhead.Rows = dimy;
    dcmhead.Columns = dimx;
    dcmhead.PixelSpacing = [pixely pixelx]';
    dcmhead.PixelAspectRatio = [1 1]';
    dcmhead.BitsAllocated = 16;
    dcmhead.BitsStored = 15;
    dcmhead.HighBit = 14;
    dcmhead.PixelRepresentation = 0;
    dcmhead.PixelPaddingValue = 0;
    dcmhead.RescaleIntercept = 0;
    dcmhead.RescaleSlope = 1;
    dcmhead.HeartRate = heart_rate;
    dcmhead.NumberOfSlices = dimz;
    dcmhead.CardiacNumberOfImages = nr_frames;
    
    
end


dicom_header = dcmhead;

end


