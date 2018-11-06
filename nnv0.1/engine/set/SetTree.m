classdef SetTree < handle
    %SetTree class for tracking reachable set
    % This SetTree is used for neural network control system
    %   Dung Tran: 11/8/2018
    
    properties
        
        S = [];
        height = 0;
        
    end
    
    methods
       % constructor
       function obj = SetTree(height)
           % @height: height of set tree
           
           if height <= 0
               error('Height of a set tree should be > 0');
           end
           
           obj.S = cell(height, 1);
           obj.height = height;
           
       end
        
       % add a reach set to some position      
       function addReachSet(obj, pos, R)
           % @pos: position of the reachable set
           % @R: an array of reachable set
           
           if pos > obj.height
               error('Posistion > limitation');
           end
           if pos <= 0
               error('Position should be >= 0');
           end
           
           obj.S{pos} = R;
        
       end
       
       % extract a reach set at some position
       function R = extractReachSet(obj, pos) 
           % @pos: position of the reachable set
           % @R: an array of reachable set
           
           if pos > obj.height
               error('Posistion > limitation');
           end
           if pos <= 0
               error('Position should be >= 0');
           end
           
           R = obj.S{pos};      
       
       end
        
       
       % get number of reachable set at a specific position
       function N = getNumOfReachSet(obj, pos)
           % @pos: position of the reachable set
           
           if pos > obj.height
               error('Posistion > limitation');
           end
           if pos <= 0
               error('Position should be >= 0');
           end
           
           N = length(obj.S{pos});           
           
       end
       
       % get total number of reachable set in the tree
       function N = getTotalNumOfReachSet(obj)
           
           N = 0;
           for i=1:obj.height
               N = N + length(obj.S{i});
           end
           
       end
       
       
       
    end
end

