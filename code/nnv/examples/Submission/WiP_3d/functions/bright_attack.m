function I = bright_attack(vol, max_pixels, threshold, noise_disturbance)
% Return a VolumeStar of a brightening attack on a few pixels

    % Initialize vars
    ct = 0; % keep track of pixels modified
    flag = 0; % determine when to stop modifying pixels
    vol = single(vol);
    at_vol = vol;

    % Create brightening attack
    for i=1:size(vol,1)
        for j=1:size(vol,2)
            for k=1:size(vol,3)
                if vol(i,j,k) < threshold
                    at_vol(i,j,k) = 255;
                    ct = ct + 1;
                    if ct >= max_pixels
                        flag = 1;
                        break;
                    end
                end
            end
            if flag == 1
                break
            end
        end
        if flag == 1
            break;
        end
    end

    % Define input set as VolumeStar
    dif_vol = -vol + at_vol;
    noise = dif_vol;
    V(:,:,:,:,1) = vol;   % center of set
    V(:,:,:,:,2) = noise; % basis vectors
    C = single([1; -1]);          % constraints
    d = single([1; -1]);          % constraints
    I = VolumeStar(V, C, d, 1-noise_disturbance, 1); % input set

end

