function dicom_header = generate_dicomheader(parameters,slice,nr_dynamics,frame,nr_frames,dimz,dimy,dimx,heart_rate,acq_dur,dcmid)

%
% GENERATES DICOM HEADER FOR EXPORT
%
% parameters = parameters from MRD file
% dynamic = current dynamic
% nr_dynamics = total number of dynamics
% frame = current frame number
% nr_frames = total number of frames
% dimz = slices / no_views2
% dimy = y dimension (phase encoding, views)
% dimx = x dimension (readout, samples)
% heart_rate = average heart rate
%
%

studyname = str2num(parameters.filename(end-9:end-6));

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
dcmhead.Format = 'DICOM';
dcmhead.FormatVersion = 3;
dcmhead.Width = dimx;
dcmhead.Height = dimy;
dcmhead.BitDepth = 15;
dcmhead.ColorType = 'grayscale';
dcmhead.FileMetaInformationGroupLength = 178;
dcmhead.FileMetaInformationVersion = uint8([0, 1])';
dcmhead.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
dcmhead.TransferSyntaxUID = '1.2.840.10008.1.2.1';
dcmhead.ImplementationClassUID = '1.2.826.0.9717382.3.0.3.6.0';
dcmhead.ImplementationVersionName = 'OFFIS_DCMTK_360';
dcmhead.SpecificCharacterSet = 'ISO_IR 100';
dcmhead.ImageType = 'ORIGINAL\PRIMARY\';
dcmhead.SOPClassUID = '1.2.840.10008.5.1.4.1.1.4';
dcmhead.StudyDate = date;
dcmhead.SeriesDate = date;
dcmhead.AcquisitionDate = date;
dcmhead.StudyTime = time;
dcmhead.SeriesTime = time;
dcmhead.AcquisitionTime = time;
dcmhead.ContentTime = time;
dcmhead.Modality = 'MR';
dcmhead.Manufacturer = 'MR Solutions Ltd';
dcmhead.InstitutionName = 'Amsterdam UMC';
dcmhead.InstitutionAddress = 'Amsterdam, Netherlands';
dcmhead.ReferringPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
dcmhead.ReferringPhysicianName.GivenName = '';
dcmhead.ReferringPhysicianName.MiddleName = '';
dcmhead.ReferringPhysicianName.NamePrefix = '';
dcmhead.ReferringPhysicianName.NameSuffix = '';
dcmhead.StationName = 'MRI Scanner';
dcmhead.StudyDescription = 'CINE';
dcmhead.SeriesDescription = '';
dcmhead.InstitutionalDepartmentName = 'Amsterdam UMC preclinical MRI';
dcmhead.PhysicianOfRecord.FamilyName = 'Amsterdam UMC preclinical MRI';
dcmhead.PhysicianOfRecord.GivenName = '';
dcmhead.PhysicianOfRecord.MiddleName = '';
dcmhead.PhysicianOfRecord.NamePrefix = '';
dcmhead.PhysicianOfRecord.NameSuffix = '';
dcmhead.PerformingPhysicianName.FamilyName = 'Amsterdam UMC preclinical MRI';
dcmhead.PerformingPhysicianName.GivenName = '';
dcmhead.PerformingPhysicianName.MiddleName = '';
dcmhead.PerformingPhysicianName.NamePrefix = '';
dcmhead.PerformingPhysicianName.NameSuffix = '';
dcmhead.PhysicianReadingStudy.FamilyName = 'AMC preclinical MRI';
dcmhead.PhysicianReadingStudy.GivenName = '';
dcmhead.PhysicianReadingStudy.MiddleName = '';
dcmhead.PhysicianReadingStudy.NamePrefix = '';
dcmhead.PhysicianReadingStudy.NameSuffix = '';
dcmhead.OperatorName.FamilyName = 'manager';
dcmhead.AdmittingDiagnosesDescription = '';
dcmhead.ManufacturerModelName = 'MRS7024';
dcmhead.ReferencedSOPClassUID = '';
dcmhead.ReferencedSOPInstanceUID = '';
dcmhead.ReferencedFrameNumber = [];  
dcmhead.DerivationDescription = '';
dcmhead.FrameType = '';
dcmhead.PatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
dcmhead.PatientID = '01';
dcmhead.PatientBirthDate = date;
dcmhead.PatientBirthTime = '';
dcmhead.PatientSex = 'F';
dcmhead.OtherPatientID = '';
dcmhead.OtherPatientName.FamilyName = 'Amsterdam UMC preclinical MRI';
dcmhead.OtherPatientName.GivenName = '';
dcmhead.OtherPatientName.MiddleName = '';
dcmhead.OtherPatientName.NamePrefix = '';
dcmhead.OtherPatientName.NameSuffix = '';
dcmhead.PatientAge = '1';
dcmhead.PatientSize = [];
dcmhead.PatientWeight = 0.0300;
dcmhead.Occupation = '';
dcmhead.AdditionalPatientHistory = '';
dcmhead.PatientComments = '';
dcmhead.BodyPartExamined = '';
dcmhead.SequenceName = 'CINE';
dcmhead.SliceThickness = parameters.SLICE_THICKNESS/dimz;
dcmhead.KVP = 0;
dcmhead.RepetitionTime = TR;     % time between frames
dcmhead.EchoTime = 2;         % approximately correct
dcmhead.InversionTime = 0;
dcmhead.NumberOfAverages = parameters.NO_AVERAGES;
dcmhead.ImagedNucleus = '1H';
dcmhead.MagneticFieldStrength = 7;
dcmhead.SpacingBetweenSlices = 0;
dcmhead.EchoTrainLength = 1;
dcmhead.DeviceSerialNumber = '0034';
dcmhead.PlateID = '';
dcmhead.SoftwareVersion = '1.0.0.0';
dcmhead.ProtocolName = 'CINE';
dcmhead.SpatialResolution = [];
dcmhead.TriggerTime = (frame-1)*TR;    % frame time = frame number times calculated TR
dcmhead.DistanceSourceToDetector = [];
dcmhead.DistanceSourceToPatient = [];
dcmhead.FieldOfViewDimensions = [parameters.FOV parameters.FOV parameters.SLICE_THICKNESS/dimz];
dcmhead.ExposureTime = [];
dcmhead.XrayTubeCurrent = [];
dcmhead.Exposure = [];
dcmhead.ExposureInuAs = [];
dcmhead.FilterType = '';
dcmhead.GeneratorPower = [];
dcmhead.CollimatorGridName = '';
dcmhead.FocalSpot = [];
dcmhead.DateOfLastCalibration = '';
dcmhead.TimeOfLastCalibration = '';
dcmhead.PlateType = '';
dcmhead.PhosphorType = '';
dcmhead.AcquisitionMatrix = uint16([dimy 0 0 dimx])';
dcmhead.FlipAngle = parameters.alpha;
dcmhead.AcquisitionDeviceProcessingDescription = '';
dcmhead.CassetteOrientation = 'PORTRAIT';
dcmhead.CassetteSize = '25CMX25CM';
dcmhead.ExposuresOnPlate = 0;
dcmhead.RelativeXrayExposure = [];
dcmhead.AcquisitionComments = '';
dcmhead.PatientPosition = 'HFS';
dcmhead.Sensitivity = [];
dcmhead.FieldOfViewOrigin = [];
dcmhead.FieldOfViewRotation = [];
dcmhead.AcquisitionDuration = acq_dur/nr_dynamics;
dcmhead.StudyInstanceUID = dcmid(1:18);
dcmhead.SeriesInstanceUID = [dcmid(1:18),'.',num2str(studyname)];
dcmhead.StudyID = '01';
dcmhead.SeriesNumber = studyname;
dcmhead.AcquisitionNumber = 1;
dcmhead.InstanceNumber = (slice-1)*nr_frames + frame;          % instance number
dcmhead.ImagePositionPatient = [(pixely/2)-parameters.FOV/2, (pixelx/2)-parameters.FOV/2, 0.0]';
dcmhead.ImageOrientationPatient = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]';
dcmhead.FrameOfReferenceUID = '';
dcmhead.TemporalPositionIdentifier = frame;
dcmhead.NumberOfTemporalPositions = nr_frames;
dcmhead.TemporalResolution = TR;
dcmhead.ImagesInAcquisition = nr_frames*dimz;
dcmhead.SliceLocation = (parameters.SLICE_THICKNESS/dimz)*((slice-1)-(round(dimz/2)));
dcmhead.ImageComments = '';
dcmhead.TemporalPositionIndex = uint32([]);
dcmhead.SamplesPerPixel = 1;
dcmhead.PhotometricInterpretation = 'MONOCHROME2';
dcmhead.PlanarConfiguration = 0;
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


dicom_header = dcmhead;


end


