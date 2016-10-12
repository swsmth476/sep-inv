classdef SepMLS
%
%  SYS: Class which contains info for a separable monotone linear system
%       problem.
%  =======================================
%  
%  
%  SYNTAX
%  ------
%     
%      SepMLS(...)
%    
%  
%  DESCRIPTION
%  -----------
%     Class to be used to search for a set of contracts for n linear, 
%     monotone subsystems. 
%  
%  INPUT
%  -----
%     
%        
%          Parameter1                  
%                       
%  
%  
%  METHODS
%  -------
%   s.trans(i, pidx, uidx)
%  
 
    properties (SetAccess=protected)
        
        sub_n; % subsystem dimensions [n1 n2 n3 n3....]'
        
        A; 
        B;
        E;
        K; % system dynamics x(k+1) = Ax(k) + Bu(k) + Ed(k)
        
        % state space partition
        xpart_res; % partition resolution (i.e. the length of a cell)
        xpart_dim; % partition dimensions (number of cells by number of cells ...)
        xpart_prod; % for partition calculations
        xpart_offset; % x-value of partition coordinate [0 0 0 0 ... 0]'
        
        % input space partition
        upart_res;
        upart_dim;
        upart_offset;
        
    end
    
    methods
        
        function ss = SepMLS(sub_n, part_res, part_dim, part_offset)
            if(size(part_dim, 2) ~= 1)
                error('Please enter partition dimensions as a column vector.');
            end
            
            if(size(sub_num, 2) ~= 1)
                error('Please enter subsystem dimensions as a column vector.');
            end
            
            if(sum(sub_num) ~= length(part_dim))
                error('Dimension of partition space should equal sum of subsystem dimensions.');
            end
            
            ss.sub_n = sub_n;
            ss.xpart_res = part_res;
            ss.xpart_dim = part_dim;
            ss.xpart_offset = part_offset;
            
            ss.xpart_prod = zeros(size(part_dim));
            for i = 1:length(part_dim)
                ss.xpart_prod(i) = prod(part_dim(1:i));
            end
        end
        
        function in = inPart(ss, pcrd)
        % note: this fn is only needed for x coordinates
            if(sum(mod(pcrd, 1)) ~= 0)
                error('Error: expected integer partition coordinates.');
            end
            if(size(pcrd) ~= size(ss.xpart_dim))
                error('Error: unexpected partition coordinate dimensions.');
            end
            in = (sum(zeros(3,1) <= pcrd <= ss.xpart_dim) == length(ss.xpart_dim));
        end
        
        function pidx = ptoi(ss, pcrd, mode)
        % partition coordinate -> partition index
            if(strcmp(mode, 'state'))
                part_prod = ss.xpart_prod;
            elseif(strcmp(mode, 'input'))
                part_prod = ss.upart_prod;
            else
                error('Error: mode specified should either be state or input.');
            end
            
            pidx = part_prod * pcrd + 1;
        end
        
        function pcrd = itop(ss, pidx, mode)
        % partition index -> partition coordinate
        
        % possibly check if partition index is out of range
        
            if(strcmp(mode, 'state'))
                part_prod = ss.xpart_prod;
            elseif(strcmp(mode, 'input'))
                part_prod = ss.upart_prod;
            else
                error('Error: mode should be a string, either "state" or "input".');
            end
        
            pcrd = zeros(size(part_prod));
            for i = 1:length(part_prod)
                pcrd(i) = floor(pidx / part_prod(i));
                pidx = mod(pidx, part_prod(i));
            end
        end
        
        function sscrd = ptox(ss, pidx, mode)
        % partition coordinate -> state space coordinate
            if(strcmp(mode, 'state'))
                part_prod = ss.xpart_prod;
            elseif(strcmp(mode, 'input'))
                part_prod = ss.upart_prod;
            else
                error('Error: mode should be a string, either "state" or "input".');
            end
            
        end
            
        function pcrd = xtop(ss, sscrd, mode)
        % state space coordinate -> partition coordinate
            
        end
        
        function next = trans(i, pidx, uidx)
            
        end
            
        
    end
    
end

