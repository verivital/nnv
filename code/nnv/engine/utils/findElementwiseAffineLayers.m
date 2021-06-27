% neelanjana pal 
function [elementwiseAffineLayers, indices] = findElementwiseAffineLayers(lgraph)
 elementwiseAffineLayers=[];
 indices = [];
     for i = length(lgraph.Layers)
         if class(lgraph.Layers(i,1)) == 'nnet.onnx.layer.ElementwiseAffineLayer'
             indices = [indices i];
             elementwiseAffineLayers = [elementwiseAffineLayers lgraph.Layers(i,1)];
         end
     end
 end