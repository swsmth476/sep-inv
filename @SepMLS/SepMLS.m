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
%      Sys('Parameter1',Value1,'Parameter2',Value2,...)
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
%          Parameter1 The name of the desired option to be     
%                     changed provided as string. The list of  
%                     available options can be obtained by     
%                     typing properties('mptopt') at the       
%                     Matlab prompt.                           
%                     Class: char                              
%          Value1     The value to be assigned to Parameter1.  
%                     Class: double or char                    
%                       
%  
%  
%  METHODS
%  -------
%   s1.transition(pidx, uidx)
%  
 
    properties (SetAccess=protected)
        
        % num subsystems
        n;
        
        % subsystem dimensions [n1 n2 n3 n3....]
        sub_n;
        
        % system dynamics x(k+1) = Ax(k) + Bu(k) + Ed(k)
        A; 
        B;
        E;
        K;
        
        % partition resolution (i.e. the length of a cell)
        part_res;
        
        % partition dimensions (number of cells by number of cells ...)
        part_dim;
        
    end
    
    methods
        
        function index = ptoi(part)
            
        end
        
        function part = xtop()
        
        function next = trans(
            
        
    end
    
end

