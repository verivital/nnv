function [images, labels, db_image_nos] = load_images(args)
% for loading MNIST and ImageNet images for ModelStar experiments (JAIR draft). Can be removed.
    
    arguments
        args.database   
        args.n               % number of images
        args.matlabnet       % if a matlab neural network is supplied, only 
            % the correctly classified images will be returned
    end
    
    database = args.database;
    n = args.n;
    matlabnet = args.matlabnet;
    
    images = {};
    labels = {};
    db_image_nos = {};
    if strcmp(database, "mnist")
        if n > 10000
            error('Maximum 10000 mnist images available.')
        end
        % Load data (no download necessary)
        digitDatasetPath = fullfile(matlabroot,'toolbox','nnet','nndemos', ...
            'nndatasets','DigitDataset');
        % Images
        imds = imageDatastore(digitDatasetPath, ...
            'IncludeSubfolders',true,'LabelSource','foldernames');
        
        numClasses = 10;
        if mod(n, numClasses) ~= 0
            error(['For MNIST, to have balanced dataset, number of images must be a multiple of ' num2str(numClasses)])
        end
        NPerClass = n/numClasses;
        
        db_img_no = 1;
        no_of_images_chosen_from_class = 0;
        class_no = 1;
        while class_no <= numClasses
            % Load one image in dataset
            [img, fileInfo] = readimage(imds, db_img_no);
            img = single(img); % change precision
            label = single(fileInfo.Label);
            
            append_image = 1;
            if ~isempty(matlabnet)
                [~, pred] = max(predict(matlabnet, img));
                if label ~= pred
                    append_image = 0;
                end
            end
            
            if append_image
                images{end + 1} = img;
                labels{end + 1} = label;
                db_image_nos{end + 1} = db_img_no;
                no_of_images_chosen_from_class = no_of_images_chosen_from_class + 1;
                if no_of_images_chosen_from_class == NPerClass
                    class_no = class_no + 1;
                    no_of_images_chosen_from_class = 0;
                    db_img_no = (class_no - 1)*1000;
                end
            end
            db_img_no = db_img_no + 1;
        end
        
    elseif strcmp(database, "imagenet")
        
        % load images from: data\examples\CNN\VGG19\weight_perturb
        folder = [nnvroot(), filesep, 'data', filesep, 'examples', filesep, 'CNN', filesep, 'VGG19', filesep, 'weight_perturb', filesep];
        files = dir(folder);
        for i = 1:length(files)
            if ~files(i).isdir
                img_orig = imread([folder files(i).name]);
                img = single(imresize(img_orig, [224 224]));
                img_class = regexp(files(i).name, '\d+', 'match');
                label = str2double(img_class{1}) + 1;
                append_image = 1;
                if ~isempty(matlabnet)
                    [~, pred] = max(predict(matlabnet, img));
                    if label ~= pred
                        append_image = 0;
                    end
                end
                
                if append_image
                    images{end + 1} = img;
                    labels{end + 1} = label;
                    if length(images) == n
                        break
                    end
                end
            end
        end
        
        if length(images) < n
            error(['Available imagenet images ' num2str(n - length(images)) ' lesser than requested.'])
        end
        
    else
        error(['Unsupported database ' database])
    end
end

