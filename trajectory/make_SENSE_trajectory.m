
kspacesize = 192;       % size of fully sampled k-space
nrcal = 30;             % size of fully sampled center of k-space
R = 2;                  % undersampling factor


nrklines = floor(kspacesize/R + nrcal - nrcal/R + 0.5);
calmin = round(kspacesize/2 - nrcal/2 + 0.5);
calmax = round(kspacesize/2 + nrcal/2 + 0.5);

k = zeros(kspacesize,1);

for i = 1:kspacesize
    
    if mod(i,R) == 0
        k(i) = i;
    else
        k(i) = NaN;
    end
    
    if i>calmin && i<calmax
        k(i) = i;
    end
    
end

k = k(~isnan(k));

% center around 0
kf = k - round(kspacesize/2 + 0.5);

% forward and back
kb = flip(kf);
kb = kb(2:end-1);
k = [kf;kb];

% length of the k-line array
disp(length(k));

fn = '/Users/Gustav/Dropbox/Reconstruction_and_analysis_tools/Matlab/MRI-apps/Retrospective5/trajectory/';

filename = strcat(fn,'sensetrajectory_L',num2str(length(k)),'_R',num2str(R),'_C',num2str(nrcal),'.txt');
fileID = fopen(filename,'w');

for i = 1:length(k)
   
    fprintf(fileID,num2str(k(i)));
    fprintf(fileID,',\n');
    
end

fclose(fileID);


%writematrix(k','/Users/Gustav/Dropbox/Reconstruction_and_analysis_tools/Matlab/MRI-apps/Retrospective5/trajectory/sensetrajectory.txt');




