classdef 2Subsys_Helper < handle
%
%  2Subsys_Helper: Helper class for 2-subsystem assumption mining.
%  =======================================
%  
%  SYNTAX
%  ------
%    
%  DESCRIPTION
%  -----------
%
%  INPUT
%  -----
%  
%  METHODS
%  -------
%    
 
    properties (SetAccess=protected)

    % disturbance partition
    dpart;
    dpart_res;
    dpart_dim;
    dpart_offset;
    
    % system interconnection
    C1;
    C2;
    E1;
    E2;
        
    end
    
    methods
       
        function 2ss = 2Subsys_Helper(dpart_res, dpart_dim, dpart_offset)
            
            if(sum(size(dpart_dim, 1) == [1 2]) ~= 2)
                error('Error: partition dimensions should be 1x2 row vector.');
            end
            
            if(sum(size(dpart_dim) == size(dpart_offset)) ~= 2)
                error('Error: partition dimension should equal offset dimension.');
            end
            
            2ss.dpart_res = dpart_res;
            2ss.dpart_dim = dpart_dim;
            2ss.dpart_offset = dpart_offset;
            
            % assumption space to search over
            dpart = -ones(dpart_dim);
            
        end
        
        function set_If(2ss, C1, C2, E1, E2)
            % maybe add more dimensions of disturbance through C matrix
            % later?
            if(size(C1,1) ~= 1 || size(C2,1) ~= 1)
                error('Only 1-dimensional disturbances supported.');
            end
            
            if(size(E1,1) ~= size(C2,2) || size(E2,1) ~= size(C1,2))
                error('Dimension mismatch between output and disturbance input.');
            end
            
            % set interconnection
            2ss.C1 = C1;
            2ss.C2 = C2;
            2ss.E1 = E1;
            2ss.E2 = E2;
        end
        
        function dcrd = to_check(
            % randomly select disturbance partition to sample
            
        end
        
        function update_search(
            % update search grid using directed spec ideas
        end
        
    end
    
end