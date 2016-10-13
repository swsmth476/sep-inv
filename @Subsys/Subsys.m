classdef Subsys < handle
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
%    for next 4 fns, mode = "state" or "input" for each partition
%    ss.overflow(ss, pcrd, mode)  check if a partition coordinate is inside the partition space
%    ss.ptoi(ss, pcrd, mode)  partition coordinate -> partition index
%    ss.itop(ss, pidx, mode)  partition index -> partition coordinate
%    mode = "upper" or "lower" for ptox, to get either upper or lower
%    corner of partition box
%    ss.ptox(ss, pidx, mode)  partition coordinate -> state space coordinate
%    ss.ptou(ss, pidx)  partition coordinate -> input space coordinate
%    ss.xtop(ss, x)  state space coordinate -> partition coordinate
%                    NOTE: only implemented for the state (not needed for input)
%    ss.setd(ss, d_set)  set disturbance assumptions for subsystem
%    ss.trans(ss, xidx, uidx)  get possible transitions
%  
 
    properties (SetAccess=protected)
        
        % subsystem dimension
        sub_n;
        
        % system dynamics x(k+1) = Ax(k) + Bu(k) + d(k)
        A;
        B;
        
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
        
        % assumption set
        d_set;
        
    end
    
    methods
        
        function ss = Subsys(sub_n, xpart_res, xpart_dim, xpart_offset, upart_res, upart_dim, upart_offset)
            
            if((size(xpart_dim, 2) ~= 1) || (size(upart_dim, 2) ~= 1))
                error('Error: partition dimensions should be a column vector.');
            end
            
            if(sum(sub_n) ~= length(xpart_dim))
                error('Error: partition dimension should equal subsystem dimension.');
            end
            
            ss.sub_n = sub_n;
            ss.xpart_res = xpart_res;
            ss.xpart_dim = xpart_dim;
            ss.xpart_offset = xpart_offset;
            
            ss.upart_res = upart_res;
            ss.upart_dim = upart_dim;
            ss.upart_offset = upart_offset;
            
            ss.xpart_prod = zeros(size(ss.xpart_dim));
            ss.xpart_prod(1) = 1;
            for i = 2:length(xpart_dim)
                ss.xpart_prod(i) = prod(xpart_dim(1:i-1));
            end
            ss.upart_prod = zeros(size(ss.upart_dim));
            ss.upart_prod(1) = 1;
            for i = 2:length(upart_dim)
                ss.upart_prod(i) = prod(upart_dim(1:i-1));
            end
        end
        
        function ovf = overflow(ss, pcrd, mode)
            
        	if(strcmp(mode, 'state'))
                part_dim = ss.xpart_dim;
                part_prod = ss.xpart_prod;
            elseif(strcmp(mode, 'input'))
                part_dim = ss.upart_dim;
                part_prod = ss.upart_prod;
            else
                error('Error: mode specified should either be "state" or "input".');
            end
            
            if(sum(mod(pcrd, 1)) ~= 0)
                error('Error: expected integer partition coordinates.');
            elseif(size(pcrd) ~= size(part_dim))
                error('Error: unexpected partition coordinate dimensions.');
            end
            
            ovf = ~((1 <= part_prod' * pcrd + 1) && (part_prod' * pcrd + 1 <= prod(part_dim)));
        end
        
        function pidx = ptoi(ss, pcrd, mode)
            if(strcmp(mode, 'state'))
                part_prod = ss.xpart_prod;
            elseif(strcmp(mode, 'input'))
                part_prod = ss.upart_prod;
            else
                error('Error: mode specified should either be "state" or "input".');
            end
            
            pidx = part_prod' * pcrd + 1;
        end
        
        function pcrd = itop(ss, pidx, mode)
            if(strcmp(mode, 'state'))
                part_prod = ss.xpart_prod;
                part_dim = ss.xpart_dim;
            elseif(strcmp(mode, 'input'))
                part_prod = ss.upart_prod;
                part_dim = ss.upart_dim;
            else
                error('Error: mode should be a string, either "state" or "input".');
            end
            
            if ~((1 <= pidx) && (pidx <= prod(part_dim)))
                error('Error: partition index out of range.');
            end
            
            pidx = pidx - 1;
            pcrd = zeros(size(part_prod));
            n = length(part_prod) + 1;
            for i = 1:length(part_prod)
                pcrd(n-i) = floor(pidx / part_prod(n-i));
                pidx = mod(pidx, part_prod(n-i));
            end
        end
        
        function x = ptox(ss, pidx, mode)
            if(strcmp(mode, 'upper'))
                x = pidx*ss.xpart_res + ss.xpart_offset + 0.5*ss.xpart_res*ones(ss.sub_n, 1);
            elseif(strcmp(mode, 'lower'))
                x = pidx*ss.xpart_res + ss.xpart_offset - 0.5*ss.xpart_res*ones(ss.sub_n, 1);
            else
                error('Error: mode should be a string, either "upper" or "lower".');
            end
        end
        
        function u = ptou(ss, pidx)
            u = pidx*ss.upart_res + ss.upart_offset;
        end
            
        function pcrd = xtop(ss, x)
            x_scaled = x./ss.xpart_res;
            pcrd = round(x_scaled);
        end
        
        function setAB(ss, A, B)
            ss.A = A;
            ss.B = B;
        end
        
        function setd(ss, d_set)
            if(size(d_set,1) ~= ss.sub_n)
               error('Error: disturbance set should have rows equal to subsystem dimension.'); 
            end
            ss.d_set = d_set;
        end
        
        function next = trans(ss, xidx, uidx)
            if(isempty(ss.d_set))
               error('Error: disturbance set has not been specified yet.'); 
            end
            
            % get state and input
            pidx = itop(ss, xidx, 'state');
            x = ptox(ss, pidx, 'upper');
            pidx = itop(ss, uidx, 'input');
            u = ptou(ss, pidx);
            
            % calculate transitions
            next = [];
            for i = 1:size(ss.d_set, 2)
                x_plus = ss.A*x + ss.B*u + ss.d_set(:, i);
                p_plus = xtop(ss, x_plus);
                if(~overflow(ss, p_plus, 'state'))
                    next = {next p_plus};
                end
            end
        end
    end
end