%to run this as a test, use results_video_to_image=runtests('test_video_to_image')
%requirements: file must start or end with test
%each test starts with two percent signs followed by the name
%shared vairables must appear before first test
%variables made by a test are not available to other tests.


%% test 1: video to image
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
