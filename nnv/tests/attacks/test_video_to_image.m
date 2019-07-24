mkdir('shuttle_images');
shuttleVideo = VideoReader('shuttle.avi');
ii = 1;

while hasFrame(shuttleVideo)
   img = readFrame(shuttleVideo);
   filename = [sprintf('%03d',ii) '.jpg'];
   fullname = fullfile('shuttle_images',filename);
   imwrite(img,fullname)    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
   ii = ii+1;
end

% play shutle video
videoFReader = vision.VideoFileReader('shuttle.avi');
videoPlayer = vision.VideoPlayer;

while ~isDone(videoFReader)
  videoFrame = videoFReader();
  videoPlayer(videoFrame);
  pause(0.1)
end