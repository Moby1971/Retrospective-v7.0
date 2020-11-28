function rpr = readrprfile(filename)


fid = fopen(filename,'r'); 
rpr = char(fread(fid,Inf,'uchar')');
fclose(fid);


end