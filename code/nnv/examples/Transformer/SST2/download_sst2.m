%% Download SST-2 Dataset
% Downloads the Stanford Sentiment Treebank (SST-2) dataset from GLUE
%
% The SST-2 dataset contains:
% - Train: ~67k sentences with binary labels (0=negative, 1=positive)
% - Dev: ~872 sentences (validation set)
% - Test: ~1.8k sentences (no labels - for leaderboard)
%
% Author: NNV Team
% Date: November 2025

function [train_data, dev_data] = download_sst2(data_dir)
    % Download and load SST-2 dataset
    %
    % Input:
    %   data_dir: Directory to store data (default: current directory)
    %
    % Output:
    %   train_data: struct with .sentences and .labels
    %   dev_data: struct with .sentences and .labels

    if nargin < 1
        data_dir = fileparts(mfilename('fullpath'));
    end

    % Create data directory
    sst2_dir = fullfile(data_dir, 'SST-2');
    if ~exist(sst2_dir, 'dir')
        mkdir(sst2_dir);
    end

    % Download URL (from GLUE benchmark)
    zip_url = 'https://dl.fbaipublicfiles.com/glue/data/SST-2.zip';
    zip_file = fullfile(data_dir, 'SST-2.zip');

    % Download if not exists
    train_file = fullfile(sst2_dir, 'train.tsv');
    dev_file = fullfile(sst2_dir, 'dev.tsv');

    if ~exist(train_file, 'file') || ~exist(dev_file, 'file')
        fprintf('Downloading SST-2 dataset...\n');

        try
            % Download zip file
            websave(zip_file, zip_url);
            fprintf('Download complete.\n');

            % Extract
            fprintf('Extracting...\n');
            unzip(zip_file, data_dir);
            fprintf('Extraction complete.\n');

            % Clean up zip
            delete(zip_file);

        catch ME
            fprintf('Download failed: %s\n', ME.message);
            fprintf('\nAlternative: Download manually from:\n');
            fprintf('  %s\n', zip_url);
            fprintf('Extract to: %s\n', sst2_dir);

            % Try to use cached/local data
            if exist(train_file, 'file')
                fprintf('Using existing local data...\n');
            else
                error('SST-2 data not available');
            end
        end
    else
        fprintf('SST-2 data already downloaded.\n');
    end

    % Load train data
    fprintf('Loading training data...\n');
    train_data = load_tsv(train_file);
    fprintf('Loaded %d training samples\n', length(train_data.sentences));

    % Load dev data
    fprintf('Loading development data...\n');
    dev_data = load_tsv(dev_file);
    fprintf('Loaded %d development samples\n', length(dev_data.sentences));
end

function data = load_tsv(filepath)
    % Load TSV file with sentence and label columns

    data = struct();
    data.sentences = {};
    data.labels = [];

    fid = fopen(filepath, 'r', 'n', 'UTF-8');
    if fid == -1
        error('Cannot open file: %s', filepath);
    end

    % Skip header
    header = fgetl(fid);

    % Read lines
    line_num = 0;
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line) && ~isempty(line)
            % Split by tab
            parts = strsplit(line, '\t');
            if length(parts) >= 2
                sentence = parts{1};
                label = str2double(parts{2});

                if ~isnan(label)
                    line_num = line_num + 1;
                    data.sentences{line_num, 1} = sentence;
                    data.labels(line_num, 1) = label;
                end
            end
        end
    end

    fclose(fid);
end
