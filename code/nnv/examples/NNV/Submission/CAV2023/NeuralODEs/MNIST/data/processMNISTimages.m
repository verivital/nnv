function X = processMNISTimages(filename)
[fileID,errmsg] = fopen(filename,'r','b');
if fileID < 0
    error(errmsg);
end

magicNum = fread(fileID,1,'int32',0,'b');
if magicNum == 2051
    fprintf('\nRead MNIST image data...\n')
end

numImages = fread(fileID,1,'int32',0,'b');
fprintf('Number of images in the dataset: %6d ...\n',numImages);
numRows = fread(fileID,1,'int32',0,'b');
numCols = fread(fileID,1,'int32',0,'b');

X = fread(fileID,inf,'unsigned char');

X = reshape(X,numCols,numRows,numImages);
X = permute(X,[2 1 3]);
% X = X./255; % Do not normalize yet
X = reshape(X, [28,28,1,size(X,3)]);
X = dlarray(X, 'SSCB');

fclose(fileID);
end

