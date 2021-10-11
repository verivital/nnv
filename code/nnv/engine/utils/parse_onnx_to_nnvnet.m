% neelanjana pal 
 function nnv_net = parse_onnx_to_nnvnet(varargin)
 %% loads an onnx model and parses it to nnv network
 %  as of now works for vnncomp2021 benchmarks:
 %                   verivital,eran,cifar2020,oval21, mnistfc 
 % need to check for others: marabou and nn4sys

 switch nargin
    case 2
      onnxfile =  varargin{1};
      ip_shape =  varargin{2};
    case 1
      onnxfile =  varargin{1};  
  end
 
  lgraph = importONNXLayers(onnxfile,'OutputLayerType','classification','ImportWeights',true);  
  [placeholderLayers, indices] = findPlaceholderLayers(lgraph);
  
  if nargin == 2
       if string(class(lgraph(1,1)))=='nnet.cnn.layer.SequenceInputLayer'
        inputlayer = lgraph(1,1);
        newInputlayer = imageInputLayer(ip_shape,'Name',inputlayer.Name,'Normalization',inputlayer.Normalization);
        lgraph(1,1)=newInputlayer;
%       elseif class(lgraph.Layers(1,1))== 'nnet.cnn.layer.SequenceInputLayer'
%         inputlayer = lgraph.Layers(1,1);  
%         newInputlayer = imageInputLayer(ip_shape,'Name',inputlayer.Name,'Normalization',inputlayer.Normalization);
%         newlayers = 
      end
  end
     
  
  if isempty(placeholderLayers) 
      if string(class(lgraph)) == 'nnet.cnn.layer.Layer'
          layers = lgraph;
      else
          layers = lgraph.Layers;
      end
     %modif_lgraph = lgraph;
     %matlabnet = assembleNetwork(lgraph);
  elseif ~isempty(placeholderLayers)
     numPH = length(indices)+1;
     idx = 1;
     for i=1:length(lgraph.Layers)
         if ~ismember(i,indices) && idx < numPH && i< indices(idx) 
            layers(i-idx+1,1) = lgraph.Layers(i,1);
         elseif ismember(i,indices) && idx < numPH && i == indices(idx)
             idx = idx + 1;
         elseif ~ismember(i,indices) && i> indices(numPH-1)
            layers(i-idx+1,1) = lgraph.Layers(i,1);
         end
     end
     
  end
  
  if nargin == 2
       if any(layers(1,1).InputSize ~= ip_shape)
        inputlayer = layers(1,1);
        newInputlayer = imageInputLayer(ip_shape,'Name',inputlayer.Name,'Normalization',inputlayer.Normalization);
        layers(1,1)=newInputlayer;
      end
  end
  for i = 1:length(layers)
     if string(layers(i,1).Name) == ''
        layers(i,1).Name = 'newNode';
     end
  end
  modif_lgraph = layerGraph(layers);
  matlabnet = assembleNetwork(modif_lgraph);
  nnv_net = CNN.parse(matlabnet);
 end
