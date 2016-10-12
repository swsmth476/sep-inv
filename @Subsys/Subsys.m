classdef Subsys
%
%  SYS: Class to represent a subsystem in a decomposable synthesis problem.
%  =======================================
%  
%  
%  SYNTAX
%  ------
%     
%      SepMLS(sub_n, xpart_res, xpart_dim, xpart_offset, upart_res, upart_dim, upart_offset)
%    
%  
%  DESCRIPTION
%  -----------
%     Class to represent a subsystem in a decomposable synthesis problem.
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
%   ss.
%  
 
    properties (SetAccess=protected)
        
        % sub_n; % subsystem dimensions [n1 n2 n3 n3....]'
        
        A; 
        B;
        E;
        K; % system dynamics x(k+1) = Ax(k) + Bu(k) + Ed(k)
        
        % state space partition
        xpart_res; % partition resolution (i.e. the length of a cell)
        xpart_dim; % partition dimensions (number of cells by number of cells ...)
        xpart_prod; % for partition calculations
        xpart_offset; % x-value of partition coordinate [0 0 0 0 ... 0]'
        % this offset corresponds to the lowest corner of our
        % partition hyperbox, i.e., it will be added to all coordinates
        
        % input space partition
        upart_res;
        upart_dim;
        upart_prod;
        upart_offset;
        
    end
    
    methods
        
        function ss = SepMLS(sub_n, xpart_res, xpart_dim, xpart_offset, upart_res, upart_dim, upart_offset)
            
            if(size(part_dim, 2) ~= 1)
                error('Error: partition dimensions should be a column vector.');
            end
            
            if(size(sub_num, 2) ~= 1)
                error('Error: subsystem dimensions should be a column vector.');
            end
            
            if(sum(sub_num) ~= length(part_dim))
                error('Error: dimension of partition space should equal sum of subsystem dimensions.');
            end
            
            ss.sub_n = sub_n;
            ss.xpart_res = xpart_res;
            ss.xpart_dim = xpart_dim;
            ss.xpart_offset = xpart_offset;
            
            ss.upart_res = upart_res;
            ss.upart_dim = upart_dim;
            ss.upart_offset = upart_offset;
            
            ss.xpart_prod = zeros(size(ss.xpart_dim));
            ss.upart_prod = zeros(size(ss.upart_dim));
            for i = 1:length(part_dim)
                ss.xpart_prod(i) = prod(part_dim(1:i));
            end
        end
        
        function ovf = overflow(ss, pcrd, mode)
        % check if a partition coordinate is inside the partition space
        	if(strcmp(mode, 'state'))
                part_dim = ss.xpart_dim;
            elseif(strcmp(mode, 'input'))
                part_dim = ss.upart_dim;
            else
                error('Error: mode specified should either be state or input.');
            end
            
            if(sum(mod(pcrd, 1)) ~= 0)
                error('Error: expected integer partition coordinates.');
            end
            if(size(pcrd) ~= size(part_dim))
                error('Error: unexpected partition coordinate dimensions.');
            end
            
            ovf = ~(sum((zeros(3,1) <= pcrd) && (pcrd <= part_dim)) == length(part_dim));
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
            if(strcmp(mode, 'state'))
                part_prod = ss.xpart_prod;
                part_dim = ss.xpart_dim;
            elseif(strcmp(mode, 'input'))
                part_prod = ss.upart_prod;
                part_dim = ss.upart_dim;
            else
                error('Error: mode should be a string, either state or input.');
            end
            
            if ~((1 <= pidx) && (pidx <= part_prod' * part_dim))
                error('Error: partition index out of range.');
            end
        
            pcrd = zeros(size(part_prod));
            for i = 1:length(part_prod)
                pcrd(i) = floor(pidx / part_prod(i));
                pidx = mod(pidx, part_prod(i));
            end
        end
        
        function sscrd = ptos(ss, pidx, mode)
        % partition coordinate -> state space coordinate
            if(strcmp(mode, 'state'))
                part_offset = ss.xpart_offset;
                part_res = ss.xpart_res;
            elseif(strcmp(mode, 'input'))
                part_offset = ss.upart_offset;
                part_res = ss.upart_res;
            else
                error('Error: mode should be a string, either state or input.');
            end
            sscrd = pidx*part_res + part_offset;
        end
            
        function pcrd = xtop(ss, x)
        % state space coordinate -> partition coordinate
        % NOTE: only implemented for the state (not needed for input)
            x_scaled = x./ss.xpart_res;
            pcrd = round(x_scaled);
        end
        
        function next = trans(ss, i, xidx, uidx)
            % get state and input
            x = ptos(ss, xidx, 'state');
            u = ptos(ss, uidx, 'input');
            
            % calculate transitions
            % ...
        end
            
        
    end
    
end

