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
    dsearch;
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
            
            % assumption space to search over
            2ss.dpart_res = dpart_res;
            2ss.dpart_dim = dpart_dim;
            2ss.dpart_offset = dpart_offset;
            
            createPart();
            dsearch = -ones(size(dpart));
            
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
        
        function createPart(2ss)
            nd = prod(2ss.dpart_dim);
            2ss.dpart = cell(1,nd);
            part = zeros(size(2ss.dpart_dim));
            2ss.dpart{1} = part;
            for i = 2:nd
                part = 2ss.incrD(part);
                2ss.dpart{i} = part;
            end
        end
        
        function dplus = incrX(ss, part_in)
            part = part_in;
            part_dim = 2ss.dpart_dim;
            n = length(2ss.dpart_dim);
            to_increment = n;
            while(part(to_increment) == part_dim(to_increment) - 1)
                to_increment = to_increment - 1;
            end
            part(to_increment) = part(to_increment) + 1;
            part(to_increment + 1:n) = 0;
            dplus = part;
        end
        
        function dcrd = find_sample()
            % randomly select disturbance partition to sample
            opt = find(dsearch == -1);
            numopt = length(opt);
            idx = randsample(numopt, 1);
            dcrd = dpart{idx};
        end
        
        function update_search(dcrd, result)
            % update search grid using directed specification ideas
            if(result == 'success')
                for i = 1:length(dsearch)
                    if(sum(dpart{i} <= dcrd) == 2)
                        dsearch(i) = 1;
                    end
                end
            else
                for i = 1:length(dsearch)
                    if(sum(dpart{i} >= dcrd) == 2)
                        dsearch(i) = 0;
                    end
                end
            end
        end
    end
end