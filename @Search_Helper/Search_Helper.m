classdef Search_Helper < handle
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
       
        function sh = Search_Helper(dpart_res, dpart_dim, dpart_offset)
            
            if(sum(size(dpart_dim, 1) == [1 2]) ~= 2)
                error('Error: partition dimensions should be 1x2 row vector.');
            end
            
            if(sum(size(dpart_dim) == size(dpart_offset)) ~= 2)
                error('Error: partition dimension should equal offset dimension.');
            end
            
            % assumption space to search over
            sh.dpart_res = dpart_res;
            sh.dpart_dim = dpart_dim;
            sh.dpart_offset = dpart_offset;
            
            createPart();
            sh.dsearch = -ones(size(dpart));
            
        end
        
        function set_If(sh, C1, C2, E1, E2)
            % maybe add more dimensions of disturbance through C matrix
            % later?
            if(size(C1,1) ~= 1 || size(C2,1) ~= 1)
                error('Only 1-dimensional disturbances supported.');
            end
            
            if(size(E1,1) ~= size(C2,2) || size(E2,1) ~= size(C1,2))
                error('Dimension mismatch between output and disturbance input.');
            end
            
            % set interconnection
            sh.C1 = C1;
            sh.C2 = C2;
            sh.E1 = E1;
            sh.E2 = E2;
            
        end
        
        function createPart(sh)
            nd = prod(sh.dpart_dim);
            sh.dpart = cell(1,nd);
            part = zeros(size(sh.dpart_dim));
            sh.dpart{1} = part;
            for i = 2:nd
                part = sh.incrD(part);
                sh.dpart{i} = part;
            end
        end
        
        function dplus = incrX(sh, part_in)
            part = part_in;
            part_dim = sh.dpart_dim;
            n = length(sh.dpart_dim);
            to_increment = n;
            while(part(to_increment) == part_dim(to_increment) - 1)
                to_increment = to_increment - 1;
            end
            part(to_increment) = part(to_increment) + 1;
            part(to_increment + 1:n) = 0;
            dplus = part;
        end
        
        function dcrd = find_sample(sh)
            % randomly select disturbance partition to sample
            opt = find(sh.dsearch == -1);
            numopt = length(opt);
            idx = randsample(numopt, 1);
            dcrd = sh.dpart{idx};
        end
        
        function update_search(sh, dcrd, result)
            % update search grid using directed specification ideas
            if(strcmp(result, 'success'))
                for i = 1:length(sh.dsearch)
                    if(sum(sh.dpart{i} <= dcrd) == 2)
                        sh.dsearch(i) = 1;
                    end
                end
            elseif(strcmp(result, 'fail'))
                for i = 1:length(sh.dsearch)
                    if(sum(sh.dpart{i} >= dcrd) == 2)
                        sh.dsearch(i) = 0;
                    end
                end
            else
                error('Invalid result. Expected "success" or "fail".');
            end
            
        end
    end
end