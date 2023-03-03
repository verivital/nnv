% load pretrained network and data
%load('M2NISTnet_92.mat');
load('M2NIST_dilatedCNN_avgpool.mat');
nnvNet = SEGNET.parse(net, 'M2NISTnet_92');
load('test_im.mat');
% attack the image by some bounded disturbance
% brightening attack for some pixels
n = size(im);
n_px_at = [5 8 11 14]; % number of pixels that are britenned
% n_px_at = [5 6]; % for testing 

% create input set
N = length(n_px_at); 
ISs(N) = ImageStar;
GrTruth = cell(1,N);

for l=1:N
    ct = 0;
    flag = 0;
    at_im = im;
    for i=1:64
        for j=1:84
            if im(i,j) > 150
                at_im(i,j) = 0;
                ct = ct + 1;
                if ct == n_px_at(l)
                    flag = 1;
                    break;
                end
            end
        end
        if flag == 1
            break;
        end
    end

    dif_im = im - at_im;

    % Perform robustness analysis
    V(:,:,:,1) = double(im);
    V(:,:,:,2) = double(dif_im);
    lb = -1;
    ub = -0.99;
    C = [1; -1];
    d = [ub; -lb];
    IS = ImageStar(V, C, d, lb, ub);
    ISs(l) = IS; 
    GrTruth{l} = im;
end

[rb, n_mis, n_rb, vr_rs] = nnvNet.verify(ISs, GrTruth, 'approx-star', 4);
nnvNet.plotSegmentationImage(im); % plot segmentation image for the input image without attack
nnvNet.plotPixelClassificationReachSet(1); % plot the pixel classification reachable set corresponding to the first input set
nnvNet.plotVerifiedOutputSet(1); % plot the first verified output set corresponding to the first input set
nnvNet.plotRobustnessStatistics(n_px_at, 'Number of Attacked Pixels'); % plot the Robustness Statistics.

