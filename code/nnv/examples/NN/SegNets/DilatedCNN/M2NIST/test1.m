
load('test_im.mat');
n_px_at = 10;
ct = 0;
flag = 0;
at_im = im;
for i=1:64
    for j=1:84
        if im(i,j) > 150
            at_im(i,j) = 0;
            ct = ct + 1;
            if ct == n_px_at
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
dif_im1 = 255*ones(64,84, 'uint8');
dif_im1 = dif_im1 - dif_im;


figure;
imshow(im);
truesize([300 300]);

figure;
imshow(dif_im1);
truesize([300 300]);

figure;
imshow(at_im);
truesize([300 300]);