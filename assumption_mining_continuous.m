% assumption mining for 2 subsystems

% subsystem dynamics
A_ii = [.1 .2; .3 .1];
B_i = [0; 1];
C_i = [1 1];
E_i = [.45; .3];

% input polyhedron
Hu = [1; -1];
hu = [.55; .65];

% subsystem 1
ls1 = LinSys(Hu, hu, A_ii, B_i, E_i, [0; 0]);

% subsystem 2
ls2 = LinSys(Hu, hu, A_ii, B_i, E_i, [0; 0]);

% search helper
sh = Search_Helper(.1, [41 41], [0 0]);
sh.set_If(E_i, E_i);

% to track results
volume_array1 = zeros(1, length(sh.dpart));
volume_array2 = zeros(1, length(sh.dpart));

for i = 1:length(sh.dpart)
    
    sample = i;

    % get assumptions
    [d1, d2, bound1, bound2] = sh.get_assumptions(sample);
    
    % set assumptions for each system
    ls1.setd([1; -1], [bound2; 0]);
    ls2.setd([1; -1], [bound1; 0]);
    
    % require that each system meets assumptions
    X1 = Polyhedron([-eye(2); C_i], [0; 0; bound1]);
    X2 = Polyhedron([-eye(2); C_i], [0; 0; bound2]);
    
    % shrink down to invariant sets
    volume1 = volume(ls1.ConInvOI(X1));
    volume2 = volume(ls2.ConInvOI(X2));
    
    % result
    if(volume1 ~= 0 && volume2 ~= 0)
        result = 'success';
    else
        result = 'fail';
    end
    
    % store results
    volume_array1(i) = volume1;
    volume_array2(i) = volume2;
    
    % update search
    sh.update_search(sample, result);
    
end