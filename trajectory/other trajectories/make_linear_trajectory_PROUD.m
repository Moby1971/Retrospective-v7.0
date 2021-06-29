







fn = '/Users/Gustav/Dropbox/Reconstruction_and_analysis_tools/Matlab/MRI-apps/Retrospective5/trajectory/';

filename = strcat(fn,'linear_32_16.txt');
fileID = fopen(filename,'w');

for k=1:32
    
    for j=1:16
   
    fprintf(fileID,num2str(k-17));
    fprintf(fileID,',');
    fprintf(fileID,num2str(j-9));
    fprintf(fileID,',\n');
    
    end
    
end

fclose(fileID);


