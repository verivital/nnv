%% Check the verification

% Define the input arguments
dsv = '';
sl = 1;
vt = '';
rt = '';
i = 1;
e = 1;

% Run the verification + create the plots
check_verify(dsv, sl, vt, rt, i, e);

%% Helper Methods
function [] = check_verify(dsVar, smpLen, vType, rType, index, eps)
    %
    % dsVar (string) : the dataset type. either "zoom_in" or "zoom_out".
    % smpLen (int)   : the length of a sample (video) in the dataset. either 4, 8, or 16.
    % vType (string) : the type of video verification. either "single_frame" or "all_frames".
    % rType (string) : the type of reachability method used. either "relax-star-area" or "approx-star".
    % index (int)    : the index into the dataset for the sample we want to analyze.
    % eps   (int)    : epsilon value for the L_inf attack.
    %

    % Check input arguments first.
    if dsVar ~= "zoom_in" && dsVar ~= "zoom_out"
        printf("dsVar argument was invalid. Must be 'zoom_in' or 'zoom_out'.")
        return
    end

    if smpLen ~= 4 && smpLen ~= 8 && smpLen ~= 16
        printf("smpLen argument was invalid. Must be 4, 8, or 16.")
        return
    end

    if vType ~= "single_frame" && vType ~= "all_frames"
        printf("vType argument was invalid. Must be 'single_frame' or 'all_frames'.")
        return
    end

    if rType ~= "relax-star-area" && rType ~= "approx-star"
        printf("rType argument was invalid. Must be 'relax-star-area' or 'approx-star'.")
        return
    end

    % Get alternative names used in file naming.
    dsVarCaps = "";
    dsVarShort = "";

    if dsVar == "zoom_in"
        dsVarCaps = "ZoomIn";
        dsVarShort = "zoomin";
    else
        dsVarCaps = "ZoomOut";
        dsVarShort = "zoomout";
    end

    fprintf("Running robustness verification on %s dataset...", dsVarCaps);

    if ~exist(sprintf("../../results/check_verify/%s/%s/", vType, dsVarCaps), "dir") % can "dir" be double quotes? it was apostrophes
        mkdir(sprintf("../../results/check_verify/%s/%s", vType, dsVarCaps));
    end

    % Load data
    data = readNPY(sprintf("../../data/%s/test/mnistvideo_%s_%df_test_data_seq.npy", dsVarCaps, dsVar, smpLen));
    labels = readNPY(sprintf("../../data/%s/test/mnistvideo_%s_test_labels_seq.npy", dsVarCaps, dsVar));

    % Preprocessing
    reshaped_data = permute(data, [1, 3, 2, 4, 5]); % to match BCSSS
    data_squeezed = squeeze(reshaped_data);
    datacopy = data_squeezed(:,:,:,:);

    % Get the sample
    sample = squeeze(datacopy(index,:,:,:));

    lb = squeeze(sample);
    ub = squeeze(sample);
    
    % Perform L_inf attack
    lb = lb - eps;
    ub = ub + eps;

    lb_min = zeros(size(sample));
    ub_max = ones(size(sample));
    lb_clip = max(lb, lb_min);
    ub_clip = min(ub, ub_max);

    VS = VolumeStar(lb_clip, ub_clip);

    % Verification
    Y_outputs = net.evaluate(sample);
    [~, yPred] = max(Y_outputs);

    % Evaluate lower and upper bounds
    LB_outputs = net.evaluate(lb_clip);
    [~, LB_Pred] = max(LB_outputs);
    UB_outputs = net.evaluate(ub_clip);
    [~, UB_Pred] = max(UB_outputs);

    % Define reachability options
    reachOptions = struct;
    reachOptions.reachMethod = rType;

    if rType == 'relax-star-area'
        reachOptions.relaxFactor = 0.5;
    end
    
    % Run verification
    t = tic;
    res_approx = net.verify_robustness(VS, reachOptions, target);
    fprintf('Robustness result with relax-star-area : %d', res_approx);

    R = net.reachSet{end};

    [lb_out, ub_out] = R.getRanges;
    lb_out = squeeze(lb_out);
    ub_out = squeeze(ub_out);

    mid_range = (lb_out + ub_out) / 2;
    range_size = ub_out - mid_range;

    x = [0 1 2 3 4 5 6 7 8 9];

    figure;
    errorbar(x, mid_range, range_size, '.');
    hold on;
    xlim([-0.5 9.5]);
    scatter(x, Y_outputs, 'x', 'MarkerEdgeColor', 'r');
    scatter(x, LB_outputs, 'x', 'MarkerEdgeColor', 'b');
    scatter(x, UB_outputs, 'x', 'MarkerEdgeColor', 'g');

    if res_approx ~= 1 && res_approx ~= 0
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star';
        res_approx = net.verify_robustness(VS, reachOptions, target);
        fprintf('Robustness result with approx-star : %d', res_approx);
    end
end
