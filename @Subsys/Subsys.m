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
%    mode = "upper" or "lower" for ptox, to get either upper or lower
%    corner of partition box
%    ss.ptox(ss, pidx, mode)  partition coordinate -> state space coordinate
%    ss.ptou(ss, pidx)  partition coordinate -> input space coordinate
%    ss.xtop(ss, x)  state space coordinate -> partition coordinate
%                    NOTE: only implemented for the state (not needed for input)
%    ss.setd(ss, d_set)  set disturbance assumptions for subsystem
%    ss.trans(ss, xidx, uidx)  get possible transitions
%    ss.makeTrans(ss)  create transition function
 
    properties (SetAccess=protected)
        
        % subsystem dimension
        sub_n;
        
        % system dynamics x(k+1) = Ax(k) + Bu(k) + d(k)
        A;
        B;
        
        % state space partition
        xpart;
        xpart_res;
        xpart_dim;
        xpart_prod;
        xpart_offset;
        
        % input space partition
        upart;
        upart_res;
        upart_dim;
        upart_prod;
        upart_offset;
        
        % assumption set
        d_set;
        
        % transition map
        tmap;
        
        % invariant set
        num_generators;
        vertices;
        inv_set;
        
    end
    
    methods
        
        function ss = Subsys(sub_n, xpart_res, xpart_dim, xpart_offset, upart_res, upart_dim, upart_offset)
            
            if((size(xpart_dim, 1) ~= 1) || (size(upart_dim, 1) ~= 1))
                error('Error: partition dimensions should be a row vector.');
            end
            
            if(sub_n ~= length(xpart_dim))
                error('Error: partition dimension should equal subsystem dimension.');
            end
            
            if(sum(size(xpart_dim) == size(xpart_offset)) ~= 2 || sum(size(upart_dim) == size(upart_offset)) ~= 2)
                error('Error: partition dimension should equal offset dimension.');
            end
            
            ss.sub_n = sub_n;
            ss.xpart_res = xpart_res;
            ss.xpart_dim = xpart_dim;
            ss.xpart_offset = xpart_offset;
            
            ss.upart_res = upart_res;
            ss.upart_dim = upart_dim;
            ss.upart_offset = upart_offset;
            
            ss.xpart_prod = zeros(size(ss.xpart_dim));
            n = length(xpart_dim);
            ss.xpart_prod(n) = 1;
            for i = 1:n-1
                ss.xpart_prod(n-i) = ss.xpart_prod(n-i+1)*ss.xpart_dim(n-i+1);
            end
            ss.upart_prod = zeros(size(ss.upart_dim));
            n = length(upart_dim);
            ss.upart_prod(n) = 1;
            for i = 1:n-1
                ss.upart_prod(n-i) = ss.upart_prod(n-i+1)*ss.upart_dim(n-i+1);
            end
            
            ss.createPart();
            ss.inv_set = zeros(size(ss.xpart));
        end
        
        function createPart(ss)
            nx = prod(ss.xpart_dim);
            ss.xpart = cell(1,nx);
            part = zeros(size(ss.xpart_dim));
            ss.xpart{1} = part;
            for i = 2:nx
                part = ss.incrX(part);
                ss.xpart{i} = part;
            end
            nu = prod(ss.upart_dim);
            ss.upart = cell(1,nu);
            part = zeros(size(ss.upart_dim));
            ss.upart{1} = part;
            for i = 2:nu
                part = ss.incrU(part);
                ss.upart{i} = part;
            end
        end
        
        function xplus = incrX(ss, part_in)
            part = part_in;
            part_dim = ss.xpart_dim;
            part_offset = ss.xpart_offset;
            n = length(ss.xpart_dim);
            to_increment = n;
            while(part(to_increment) == part_dim(to_increment) - 1)
                to_increment = to_increment - 1;
            end
            part(to_increment) = part(to_increment) + 1;
            part(to_increment + 1:n) = 0;
            xplus = part;
        end
        
        function uplus = incrU(ss, part_in)
            part = part_in;
            part_dim = ss.upart_dim;
            part_offset = ss.upart_offset;
            n = length(ss.upart_dim);
            to_increment = n;
            while(part(to_increment) == part_dim(to_increment) - 1)
                to_increment = to_increment - 1;
            end
            part(to_increment) = part(to_increment) + 1;
            part(to_increment + 1:n) = 0;
            uplus = part;
        end
        
        function ovf = overflow(ss, pcrd, mode)
        	if(strcmp(mode, 'state')) %#ok<ALIGN>
                part_dim = ss.xpart_dim;
                part_offset = ss.xpart_offset;
            elseif(strcmp(mode, 'input'))
                part_dim = ss.upart_dim;
                part_offset = ss.upart_offset;
            else
                error('Error: mode specified should either be "state" or "input".');
            end
            
            if(sum(mod(pcrd, 1)) ~= 0)
                error('Error: expected integer partition coordinates.');
            elseif(size(pcrd) ~= size(part_dim))
                error('Error: unexpected partition coordinate dimensions.');
            end
            
            n = length(part_dim);
            ovf = (sum(zeros(size(part_dim)) <= pcrd) ~= n) ...
                            || (sum(pcrd < part_dim) ~= n);
        end
        
        function pidx = ptoi(ss, pcrd, mode)
            if(strcmp(mode, 'state'))
                if(ss.overflow(pcrd, 'state'))
                    error('Error: coordinate is out of range of state partition.');
                end
                part_prod = ss.xpart_prod;
            elseif(strcmp(mode, 'input'))
                if(ss.overflow(pcrd, 'input'))
                    error('Error: coordinate is out of range of input partition.');
                end
                part_prod = ss.upart_prod;
            else
                error('Error: mode specified should either be "state" or "input".');
            end
            
            if(sum(size(part_prod) == size(pcrd)) ~= 2)
                error('Error: unexpected coordinate dimension.');
            end
            
            pidx = part_prod * pcrd' + 1;
        end
        
        function x = ptox(ss, pidx, mode)
            if(ss.overflow(pidx, 'state'))
                error('Error: coordinate is out of state partition.');
            end
            
            if(strcmp(mode, 'upper'))
                x = pidx.*ss.xpart_res + ss.xpart_offset + 0.5.*ss.xpart_res.*ones(1, ss.sub_n);
            elseif(strcmp(mode, 'lower'))
                x = pidx.*ss.xpart_res + ss.xpart_offset - 0.5.*ss.xpart_res.*ones(1, ss.sub_n);
            else
                error('Error: mode should be a string, either "upper" or "lower".');
            end
        end
        
        function u = ptou(ss, pidx)
            if(ss.overflow(pidx, 'input'))
                error('Error: coordinate is out of input partition.');
            end
            u = pidx.*ss.upart_res + ss.upart_offset;
        end
            
        function pcrd = xtop(ss, x)
            x_scaled = x./ss.xpart_res;
            % to ensure that we round DOWN from half-coordinates
            % (according to Sam's transition properties)
            
            % for floating point errors
            threshold = ss.xpart_res*1e-4;
            
            x_scaled(abs(mod(x_scaled, 1) - 0.5) < threshold) = ...
                x_scaled(abs(mod(x_scaled, 1) - 0.5) < threshold) - 0.5;
            pcrd = round(x_scaled);
        end
        
        function setAB(ss, A, B)
            if(sum(size(A) == [ss.sub_n ss.sub_n]) ~= 2)
                error(['wrong dimensions of A, found ', num2str(size(A)), ' expected ', num2str([ss.sub_n ss.sub_n])]);
            end
            if((size(B,1) ~= ss.sub_n) || (size(B,2) ~= size(ss.upart_dim,2)))
                error(['wrong dimensions of B, found ', num2str(size(B)), ' expected ', num2str([ss.sub_n size(ss.upart_dim,2)])]);
            end
            ss.A = A;
            ss.B = B;
        end
        
        function setd(ss, d_set)
            % add all upper bounds of disturbance as well as lower bounds
            % of disturbance
            if(size(d_set,2) ~= ss.sub_n)
               error('Error: disturbance set should have cols equal to subsystem dimension.'); 
            end
            ss.d_set = d_set;
        end
        
        function next = trans(ss, xidx, uidx)
            if(isempty(ss.d_set))
               error('Error: disturbance set has not been specified yet.'); 
            end
            
            if(isempty(ss.A) || isempty(ss.B))
               error('Error: dynamics not specified.');
            end
            
            pidx = ss.xpart{xidx};
            x = ptox(ss, pidx, 'upper');
            pidx = ss.upart{uidx};
            u = ptou(ss, pidx);
            
            next = [];
            for i = 1:size(ss.d_set, 1)
                x_plus = (ss.A*x' + ss.B*u' + ss.d_set(i, :)')';
                p_plus = xtop(ss, x_plus);
                if(~overflow(ss, p_plus, 'state'))
                    next{i} = ptoi(ss, p_plus, 'state');
                else
                    % it overflowed, add indicator of this
                    next{i} = -1;
                end
            end
        end
        
        function getTrans(ss)
            ss.tmap = cell(length(ss.xpart), length(ss.upart));
            for i = 1:length(ss.xpart)
                for j = 1:length(ss.upart)
                    ss.tmap{i, j} = ss.trans(i, j);
                end
            end
        end
        
        function setGen(ss, gen)
            % maybe check to make sure gen makes geometric sense?
            ss.num_generators = gen;
        end
        
        function setInv(ss, vert_in)
            if(isempty(ss.num_generators))
                error('Error: number of invariant set generators unset.');
            end
            if(sum(size(vert_in) == [ss.num_generators ss.sub_n]) ~= 2)
                error(['wrong dimensions of input, found ', num2str(size(vert_in)), ' expected ', num2str([ss.num_generators ss.sub_n])]);
            end
            ss.vertices = vert_in;
            ss.updateInv();
        end
        
        function in = inside(ss, xidx)
           in = 0;
           for i = 1:ss.num_generators
               if(sum(ss.xpart{xidx} <= ss.vertices(i,:)) == ss.sub_n)
                   in = 1;
                   break;
               end
           end
        end
        
        function updateInv(ss)
            for i = 1:length(ss.xpart)
                ss.inv_set(i) = ss.inside(i);
            end
        end
        
        function forall = isSafe(ss, x)
            % x, u, and d are all indices
            % ss.inside(ss.tmap{x, u}{d})
            % \forall d, \exists u s.t. above holds
            forall = 1;
            for d = 1:size(ss.d_set, 1)
                exists = 0;
                for u = 1:length(ss.upart)
                    if(ss.tmap{x, u}{d} ~= -1 && ss.inv_set(ss.tmap{x, u}{d}))
                        exists = 1;
                        % disp(['Transition from ', num2str(x), ' to ', ...
                            % num2str(ss.tmap{x, u}{d}), ' via ', num2str(u)]);
                        break;
                    end
                end
                if(~exists)
                   forall = 0;
                   break;
                end
            end
        end
        
        function forall = verifyInv(ss)
            % check every state in the invariant set
            forall = 1;
            for i = find(ss.inv_set)
                if(~ss.isSafe(i))
                    forall = 0;
                    break;
                end
            end
        end
    end
end