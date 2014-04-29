[file,path]=uigetfile('*');

vid = VideoReader([path,file]);
vido = VideoWriter([path,'slowed']);
vido.FrameRate = 2;
open(vido);
for i = 1:vid.NumberOfFrames
    frame = read(vid,i);
    writeVideo(vido,frame);
end
close(vido)