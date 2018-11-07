classdef SetTree < handle
    %SetTree class for tracking reachable set
    % This SetTree is used for neural network control system
    %   Dung Tran: 11/8/2018
    
    properties
        
        S = []; % reach set tree
        fb = []; % feedback set tree
        height = 0; % height of the set tree
        
    end
    
    methods
       % constructor
       function obj = SetTree(height)
           % @height: height of set tree
           
           if height <= 0
               error('Height of a set tree should be > 0');
           end
           
           obj.S = cell(height, 1);
           obj.fb = cell(height, 1);
           obj.height = height;
           
       end
        
       % add a reach set to some position      
       function addReachSet(obj, R, pos)
           % @pos: position of the reachable set
           % @R: an array of reachable set
           
           if pos > obj.height
               error('Posistion > limitation');
           end
           if pos <= 0
               error('Position should be > 0');
           end          
           
           obj.S{pos} = R; % update set tree
           
           % update feedback set tree
           if pos == 1
               obj.fb{pos} = cell(1,1);
               obj.fb{pos}{1, 1} = R;
           elseif pos > 1
               l = length(obj.fb{pos - 1});
               p = length(R);
               obj.fb{pos} = cell(1, l*p); 
               
               for k=1:l*p
                   i = fix((k-1)/p) + 1;
                   j = k - (i - 1) * p;
                   old_fb = obj.fb{pos-1}{1, i};                   
                   obj.fb{pos}{1, k} = [old_fb R(j)]; 
               end
               
           end
        
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
       
       % extract feedback reach set
       function fb_R = extract_fb_ReachSet(obj, pos)
           if pos > obj.height
               error('Posistion > limitation');
           end
           if pos <= 0
               error('Position should be >= 0');
           end
           
           fb_R = obj.fb{pos}; % a cell of cell   
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
       
       % flattening SetTree
       function R = flatten(obj)
           
           R = [];
           for i=1:obj.height
               R = [R obj.S{i}];
           end
           
       end
            
    end
end

