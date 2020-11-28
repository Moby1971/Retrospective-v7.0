function footer = readmrdfooter(mrdfile)

% Read the footer of an MRD file



% Read information from the header and footer first
fid = fopen(mrdfile,'r');  
val = fread(fid,4,'int32');
xdim = val(1); 
ydim = val(2); 
zdim = val(3);
dim4 = val(4);
fseek(fid,18,'bof');
datatype=fread(fid,1, 'uint16');
datatype = dec2hex(datatype);
fseek(fid,152,'bof');
val = fread(fid,2, 'int32');
dim5 = val(1);
dim6 = val(2);
no_samples = xdim;  
no_views = ydim;    
no_views_2 = zdim;  
no_slices = dim4;
no_echoes = dim5;
no_expts = dim6;


% Determine datatype
if size(datatype,2)>1
    onlydatatype = datatype(2);
    iscomplex = 2;
else
    onlydatatype = datatype(1);
    iscomplex = 1;
end
switch onlydatatype
    case '0' 
        dataformat = 'uchar';   datasize = 1; % size in bytes
    case '1' 
        dataformat = 'schar';   datasize = 1; % size in bytes
    case '2' 
        dataformat = 'short';   datasize = 2; % size in bytes
    case '3' 
        dataformat = 'int16';   datasize = 2; % size in bytes
    case '4' 
        dataformat = 'int32';   datasize = 4; % size in bytes
    case '5' 
        dataformat = 'float32'; datasize = 4; % size in bytes
    case '6' 
        dataformat = 'double';  datasize = 8; % size in bytes
    otherwise
        dataformat = 'int32';   datasize = 4; % size in bytes
end


% Fast forward to beginning of footer
fseek(fid,512,'bof');
num2read = no_expts*no_echoes*no_slices*no_views_2*no_views*no_samples*iscomplex;
fseek(fid,num2read*datasize,'cof');


% Read the footer
% sample_filename = char(fread(fid,120,'uchar')');
footer = char(fread(fid,Inf,'uchar')');
fclose(fid);


end