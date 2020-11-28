function export_gif(gifexportpath,movie,tag,window,level)
% exports movie to animated gif


% dimensions
[number_of_frames,~,~,number_of_slices,number_of_dynamics] = size(movie);

% create folder if not exist, and clear
folder_name = [gifexportpath,[filesep,'CINE_GIFS_',num2str(number_of_frames),'_',num2str(number_of_slices),'_',num2str(number_of_dynamics),'_',tag]];
if (~exist(folder_name, 'dir')); mkdir(folder_name); end
delete([folder_name,filesep,'*']);

% one heartbeat is 1 second in the movie
delay_time = 1/number_of_frames;  

% scale from 0 to 255
window = window*255/max(movie(:));
level = level*255/max(movie(:));
movie = movie*255/max(movie(:));

% window and level
movie = (255/window)*(movie - level + window/2);
movie(movie < 0) = 0;
movie(movie > 255) = 255;


for i = 1:number_of_slices
    
    slice = ['0',num2str(i)];
    slice = slice(end-1:end);
    
    for j = 1:number_of_dynamics
        
        dynamic = ['00',num2str(j)];
        dynamic = dynamic(end-2:end);
        
        for idx = 1:number_of_frames
            
            if idx == 1
                imwrite(uint8(squeeze(movie(idx,:,:,i,j))),[folder_name,filesep,'movie_',tag,'_slice',slice,'_dynamic',dynamic,'.gif'],'DelayTime',delay_time,'LoopCount',inf);
            else
                imwrite(uint8(squeeze(movie(idx,:,:,i,j))),[folder_name,filesep,'movie_',tag,'_slice',slice,'_dynamic',dynamic,'.gif'],'WriteMode','append','DelayTime',delay_time);
            end
            
        end
        
    end
    
end

end                 