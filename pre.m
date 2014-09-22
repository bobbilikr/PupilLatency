
vid=videoinput('winvideo', 1,'YUY2_320x240');          % Video Parameters
set(vid,'ReturnedColorSpace','grayscale');      % acquire in greyscale
figure;
hold on;
while(1)
    imshow(getsnapshot(vid));
    pause(0.01);
end% Preview