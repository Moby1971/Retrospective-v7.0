function [rawdata,parameters] = importb(import_path)


% test data
%import_path = '/Users/gustav/Dropbox/Reconstruction_and_analysis_tools/Matlab/MRI-apps example data/Retrospective/Bruker data/VolumeResonator_PV_6.0.1/short_axis/22_short_scan/';
%import_path = '/Users/gustav/Dropbox/Reconstruction_and_analysis_tools/Matlab/MRI-apps example data/Retrospective/Bruker data/4channel_cardio_PV_6.0.1/short_axis/22_long_scan/';
%import_path = '/Users/gustav/Dropbox/Reconstruction_and_analysis_tools/Matlab/MRI-apps example data/Retrospective/Bruker data/VolumeResonator_PV_5.0/short_axis_sperate_slice_nav/7/';


% Parameters
info1 = jcampread(strcat(import_path,'acqp'));
info2 = jcampread(strcat(import_path,'method'));



% Slices
parameters.NO_SLICES = str2num(info1.NSLICES);
parameters.SLICE_THICKNESS = str2num(info2.pvm.slicethick) * parameters.NO_SLICES;



% Matrix in readout direction
parameters.NO_SAMPLES = info1.acq.size(1)/2;
if isfield(info2.pvm,"matrix")
    parameters.NO_VIEWS = info2.pvm.encmatrix(1);
end



% Matrix in phase encoding direction
parameters.NO_VIEWS = info1.acq.size(2);
parameters.PHASE_ORIENTATION = 1;
if isfield(info2.pvm,"matrix")
    parameters.NO_VIEWS = info2.pvm.encmatrix(2);
end



% Matrix in 2nd phase encoding direction
parameters.NO_VIEWS_2 = 1;
parameters.pe2_centric_on = 0;



% FOV
parameters.FOV = info1.acq.fov(1)*10;



% Sequence parameters
parameters.tr = info1.acq.repetition_time;
parameters.te = info1.acq.echo_time;
parameters.alpha = str2num(info1.acq.flip_angle);
parameters.NO_ECHOES = 1;
parameters.NO_AVERAGES = str2num(info1.NA);



% Other parameters
parameters.date = datetime;
parameters.nucleus = '1H';
parameters.PPL = 'Navigator Sequence';
parameters.filename = 'Proton';
parameters.field_strength = str2num(info1.BF1)/42.58;
parameters.filename = 111;
parameters.pe1_order = 2;
parameters.radial_on = 0;
parameters.slice_nav = 1;



% Number of k-space repetitions
if str2num(info1.NR) > 1
    parameters.EXPERIMENT_ARRAY = str2num(info1.NR);
else
    parameters.EXPERIMENT_ARRAY = info1.acq.size(2)/info2.pvm.encmatrix(2);
end



% Number of dynamics (not sure this is needed)
if isfield(info2,"Concatenations")
    parameters.EXPERIMENT_ARRAY = parameters.EXPERIMENT_ARRAY * str2num(info2.Concatenations);
end



% Number of navigator points
if isfield(info2.pvm,"navpoints")
    parameters.no_samples_nav = str2num(info2.pvm.navpoints);
else
    parameters.no_samples_nav = str2num(info2.NavSize)/2;
end



% Number of receiver coils
parameters.no_coils = str2num(info2.pvm.encnreceivers);



% Trajectory 1st phase encoding direction
if isfield(info2.pvm,'ppggradamparray1')
    parameters.gp_var_mul = round(-info2.pvm.ppggradamparray1 * (parameters.NO_VIEWS/2-0.5));
    parameters.pe1_order = 3;
elseif isfield(info2.pvm,'encvalues1')
    parameters.gp_var_mul = round(-info2.pvm.encvalues1 * (parameters.NO_VIEWS/2-0.5));
    parameters.pe1_order = 3;
else
    % assume zigzag
    parameters.pe1_order = 2;
end



% Data type
datatype = 'int32';
if isfield(info1.acq,'word_size')
    if strcmp(info1.acq.word_size,'_32_BIT') 
        datatype = 'int32';
    end
     if strcmp(info1.acq.word_size,'_16_BIT') 
        datatype = 'int16';
    end
end



% Read data
if isfile(strcat(import_path,'fid.orig'))
    fileID = fopen(strcat(import_path,'fid.orig'));
else
    fileID = fopen(strcat(import_path,'rawdata.job0'));
end
data = fread(fileID,datatype);
fclose(fileID);
kreal = data(1:2:end);
kim = data(2:2:end);
kspace = kreal + 1j*kim;



% Read navigator
if isfile(strcat(import_path,'fid.NavFid'))
    fileID = fopen(strcat(import_path,'fid.NavFid'));
else
    fileID = fopen(strcat(import_path,'rawdata.job1'));
end
navdata = fread(fileID,datatype);
fclose(fileID);
kreal = navdata(1:2:end);
kim = navdata(2:2:end);
navkspace = kreal + 1j*kim;


% Phase offset
if isfield(info1.acq,'phase1_offset')
    parameters.pixelshift1 = round(-parameters.NO_VIEWS*info1.acq.phase1_offset/parameters.FOV);
end


% 2D data
if strcmp(info2.pvm.spatdimenum,"2D") || strcmp(info2.pvm.spatdimenum,"<2D>")
    
    kspace = reshape(kspace,parameters.NO_SLICES,parameters.NO_SAMPLES,parameters.no_coils,parameters.NO_VIEWS,parameters.EXPERIMENT_ARRAY);
    kspace = permute(kspace,[3,5,1,4,2]); % nc, nr, ns, np, nf
    
    navkspace = reshape(navkspace,parameters.NO_SLICES,parameters.no_samples_nav,parameters.no_coils,parameters.NO_VIEWS,parameters.EXPERIMENT_ARRAY);
    navkspace = permute(navkspace,[3,5,1,4,2]);
    
    kspacer = zeros(parameters.no_coils,parameters.EXPERIMENT_ARRAY,parameters.NO_SLICES,parameters.NO_VIEWS,34);
    
    raw = cat(5,navkspace,kspacer,kspace);
    for i = 1:parameters.no_coils
        rawdata{i} = squeeze(raw(i,:,:,:,:));
        rawdata{i} = reshape(rawdata{i},parameters.EXPERIMENT_ARRAY,parameters.NO_SLICES,parameters.NO_VIEWS,parameters.NO_SAMPLES+34+parameters.no_samples_nav);
    end
    
end



% 3D data (NOT EXTENSIVELY TESTED)
if strcmp(info2.pvm.spatdimenum,"3D") || strcmp(info2.pvm.spatdimenum,"<3D>")
    
    parameters.NO_VIEWS_2 = info1.acq.size(3);
    if isfield(info2.pvm,"matrix")
        parameters.NO_VIEWS = info2.pvm.encmatrix(3);
    end
    
    % Phase offset
    if isfield(info1.acq,'phase2_offset')
        parameters.pixelshift2 = round(-parameters.NO_VIEWS_2*info1.acq.phase2_offset/parameters.FOV);
    end
    
    parameters.SLICE_THICKNESS = str2num(info2.pvm.slicethick);
    parameters.pe2_centric_on = 2;
    
    if isfield(info2.pvm,"encsteps2")
        parameters.pe2_traj = info2.pvm.encsteps2;
    end
    
    if isfield(info2.pvm,'encvalues2')
        parameters.pe2_traj = round(info2.pvm.encvalues2 * (parameters.NO_VIEWS_2/2-0.5));
    end
    
    kspace = reshape(kspace,parameters.no_coils,parameters.NO_SAMPLES,parameters.NO_VIEWS,parameters.NO_VIEWS_2,parameters.EXPERIMENT_ARRAY);
    kspace = permute(kspace,[1,5,4,3,2]);
    
    navkspace = reshape(navkspace,parameters.no_coils,parameters.no_samples_nav,parameters.NO_VIEWS,parameters.NO_VIEWS_2,parameters.EXPERIMENT_ARRAY);
    navkspace = permute(navkspace,[1,5,4,3,2]);
    
    kspacer = zeros(parameters.no_coils,parameters.EXPERIMENT_ARRAY,parameters.NO_VIEWS_2,parameters.NO_VIEWS,34);
    
    raw = cat(5,navkspace,kspacer,kspace);
    for i = 1:parameters.no_coils
        rawdata{i} = squeeze(raw(i,:,:,:,:));
        rawdata{i} = reshape(rawdata{i},parameters.EXPERIMENT_ARRAY,parameters.NO_VIEWS_2,parameters.NO_VIEWS,parameters.NO_SAMPLES+34+parameters.no_samples_nav);
    end
    
end





