imaqhwinfo                              % Stands for image acquisition hardware info
cam=imaqhwinfo;
cam.InstalledAdaptors                   % Shows the data about the Installed adaptors
dev_info = imaqhwinfo('winvideo',1)     % info about the Device
vid=videoinput('winvideo', 1,'YUY2_320x240');          % Video Parameters
preview(vid);                           % Preview