
vid=videoinput('winvideo', 1,'YUY2_320x240');          % Video Parameters
set(vid,'ReturnedColorSpace','grayscale');      % acquire in greyscale
triggerconfig(vid, 'manual');
start(vid);
gcf=figure;
set(gcf,'CloseRequestFcn',@my_closefcn)
hold on;
closeflag=1;
while(closeflag)
    imshow(getsnapshot(vid));
    pause(0.01);
end% Preview

function my_closefcn()
    closeflag=0;
    delete(gcf);
end