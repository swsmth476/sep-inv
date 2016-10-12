classdef Subsys
%
%  SUBSYS: Class to represent a subsystem in a decomposable synthesis problem.
%  =======================================
%  
%  SYNTAX
%  ------
%    Subsys(sub_n, xpart_res, xpart_dim, xpart_offset, upart_res, upart_dim, upart_offset)
%    
%  DESCRIPTION
%  -----------
%    Contains a state space partition for both the state and input of the
%    subsystem being considered. Can be used to compute a transition 
%    function for the subsystem under consideration.
%
%  INPUT
%  -----
%    xpart_res      Resolution of the state space partition (i.e. the length of a cell)
%    xpart_dim      Dimension of the state space partition (number of cells by number of cells ...)
%    xpart_offset   Corresponds to the lowest corner of our partition hyperbox
%                   with partition coordinate [0 0 0 ... 0]. Will be added to all
%                   coordinates
%  
%    upart_...      Similarly defined as above for the input space
%                   partition
%  
%  METHODS
%  -------
%   ss.overflow(ss, pcrd, mode)  check if a partition coordinate is inside the partition space
%   ss.ptoi(ss, pcrd, mode)  partition coordinate -> partition index
%   ss.itop(ss, pidx, mode)  partition index -> partition coordinate
%   ss.ptos(ss, pidx, mode)  partition coordinate -> state space coordinate
%   ss.xtop(ss, x)  state space coordinate -> partition coordinate
%                          NOTE: only implemented for the state (not needed for input)
%  
 
    properties (SetAccess=protected)
        
        % system dynamics x(k+1) = Ax(k) + Bu(k) + Ed(k)
        A; 
        B;
        E;
        K;
        
        % state space partition
        xpart_res;
        xpart_dim;
        xpart_prod;
        xpart_offset;
        
        % input space partition
        upart_res;
        upart_dim;
        upart_prod;
        upart_offset;
        
    end
    
    methods
        
        function ss = SepMLS(xpart_res, xpart_dim, xpart_offset, upart_res, upart_dim, upart_offset)
            
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
            x_scaled = x./ss.xpart_res;
            pcrd = round(x_scaled);
        end
        
        function next = trans(ss, xidx, uidx)
            % get state and input
            x = ptos(ss, xidx, 'state');
            u = ptos(ss, uidx, 'input');
            
            % calculate transitions
            % ...
        end
            
        
    end
    
end

