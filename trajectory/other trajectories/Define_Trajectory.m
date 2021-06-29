close all
clear all
clc

%% Get data
% [dataraw,dimensions,parameters] = Get_mrd_3D4('/home/cylandegent/lood_storage/divh/BMEPH/Gustav/Users/cylandegent/05. Radial data 06-05-2020/11710/11710_0.MRD','seq','cen'); % phantom 06-05-2020
% 
% % rearrange data to y,x
% if ndims(dataraw) == 3
%     dataraw = permute(dataraw,[3 2 1]);
% else
%     dataraw = permute(dataraw,[2 1]);
% end
% 
% [dimx,dimy,ns] = size(dataraw);

dimx = 256;
dimy = 128;
ns = 1;

%% Parameters
radial = 360;
U_spokes = 1;

movie = 1;
half_spoke = 0;

%% Trajectory
traj = zeros(3,dimx,dimy);
for n = 1:dimy
    
    % caluculate angle
    if radial == 180
        a(n) = (n-1)*180/dimy;
    elseif ~ U_spokes
        a(n) = (n-1)*360/dimy;
    elseif (n-1) < dimy/2
        a(n) = (n-1)*360/dimy;
    else
        a(n) = (n-1)*360/dimy + 180/dimy;
    end
    
    % calculate x and y values
    for m = 1:dimx
        r(m,n) = (dimx-1)/2 -(m-1);
        x = r(m,n)*-cos(a(n)*(pi/180));
        y = r(m,n)*sin(a(n)*(pi/180));
        
        traj(1,m,n) = x;
        traj(2,m,n) = y;
    end
end

%% Make a movie or plot from data
if movie
    
    %initializing movie
    outfile = sprintf('%s','movie');
    mov = VideoWriter('movie.avi') %,'compression','Cinepak');
    mov.FrameRate = 10;
    mov.Quality = 100;
    
    open(mov)
    
    %initializing figures
    fig2 = figure(2);
    h = fig2;
    winsize = get(fig2,'Position');
    winsize(1:2) = [0 0];
    title('simulation');
    
    ColorMap = jet(dimx);
    kmax = (dimx-1)/2 +15;
    for n = 1:dimy
        
        if half_spoke
            plot(traj(1,1:end/2,n),traj(2,1:end/2,n),'color',ColorMap(n*2,:));
        else
            if n<=dimy/2
                plot(traj(1,:,n),traj(2,:,n),'color','r');
            else
                plot(traj(1,:,n),traj(2,:,n),'color',ColorMap(5,:));
            end
        end
        
        axis([-kmax kmax -kmax kmax]);
        drawnow
        hold on
        h = figure(2);
        F = getframe(h);
        writeVideo(mov,F);
    end
    
    close(mov);
    disp('FINISHED MOVIE');
    
else
    figure()
    plot((squeeze(traj(1,:,:))), (squeeze(traj(2,:,:))))
    title('radial k-space trajectory')
end
