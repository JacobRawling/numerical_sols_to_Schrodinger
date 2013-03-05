%% Open avi file
writerObj = VideoWriter('a wavefunc');
open(writerObj);
G = newplot;

%% Add one frame at a time
for x = 0:10:14999%start:increments:end
    format LONGG;
    num = num2str(x,5);
    name = strcat('time',num);
    name = strcat(name,'.txt');
    data = load(name);
    
    plot(data(1:400,1:1),data(1:400,2:2),data(1:400,1:1),data(1:400,3:3));
    axis([0 400 0 1]);
    G = getframe;
    drawnow 
    % Fetch frame
    writeVideo(writerObj,G);
end

%% Close video
close(writerObj);
