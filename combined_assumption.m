% assumption mining for 2 subsystems

% subsystem dynamics
A_ii = [.1 .2; .3 .1];
B_i = [0; 1];
C_i = [1 1];
E_i = [.45; .3];

% subsystem 1
ss1 = Subsys(2, .1, [21 21], [0 0], .1, 13, -.65);
ss1.setAB(A_ii, B_i);

% subsystem 2
ss2 = Subsys(2, .1, [21 21], [0 0], .1, 13, -.65);
ss2.setAB(A_ii, B_i);

% search helper
sh = Search_Helper(.1, [41 41], [0 0]);
sh.set_If(E_i, E_i);

% input polyhedron
Hu = [1; -1];
hu = [.55; .65];

% subsystem 1
ls1 = LinSys(Hu, hu, A_ii, B_i, E_i, [0; 0]);

% subsystem 2
ls2 = LinSys(Hu, hu, A_ii, B_i, E_i, [0; 0]);

% to track results
discrete_array1 = zeros(1, length(sh.dpart));
discrete_array2 = zeros(1, length(sh.dpart));
discrete_array3 = zeros(1, length(sh.dpart));
continuous_array1 = zeros(1, length(sh.dpart));
continuous_array2 = zeros(1, length(sh.dpart));
continuous_array3 = zeros(1, length(sh.dpart));
cWin = zeros(1, length(sh.dpart));
dWin = zeros(1, length(sh.dpart));
tie = zeros(1, length(sh.dpart));

for i = 1:length(sh.dpart)
    
    sample = i

    % get assumptions
    [d1, d2, bound1, bound2] = sh.get_assumptions(sample);
    
    % set assumptions for each system
    ss1.setd(d2);
    ss2.setd(d1);
    ss1.getTrans();
    ss2.getTrans();
    
    % require that each system meets assumptions
    ss1.setSafe(C_i, bound1);
    ss2.setSafe(C_i, bound2);
    ss1.updateInv();
    ss2.updateInv();
    
    % shrink down to invariant sets
    volume1 = ss1.ConInvOI();
    volume2 = ss2.ConInvOI();
    
    % store results
    discrete_array1(i) = volume1;
    discrete_array2(i) = volume2;
    
    if(volume1 == 0 || volume2 == 0)
        discrete_array3(i) = 0;
    else
        discrete_array3(i) = volume1 + volume2;
    end
    
    if(volume1 ~= 0 && volume2 ~= 0)
        discrete_result = 1;
    else
        discrete_result = 0;
    end
    
    % set assumptions for each system
    ls1.setd([1; -1], [bound2; 0]);
    ls2.setd([1; -1], [bound1; 0]);
    
    % require that each system meets assumptions
    X1 = Polyhedron([-eye(2); C_i], [0; 0; bound1]);
    X2 = Polyhedron([-eye(2); C_i], [0; 0; bound2]);
    
    % shrink down to invariant sets
    [C1, iter1] = ls1.ConInvOI(X1);
    [C2, iter2] = ls2.ConInvOI(X2);
    
    volume1 = volume(C1);
    volume2 = volume(C2);
    
    continuous_array1(i) = volume1;
    continuous_array2(i) = volume2;
    
    if(volume1 == 0 || volume2 == 0)
        continuous_array3(i) = 0;
    else
        continuous_array3(i) = volume1 + volume2;
    end
    
    if(volume1 ~= 0 && volume2 ~= 0)
        continuous_result = 1;
    else
        continuous_result = 0;
    end
    
    if(discrete_result == 1 && continuous_result == 0)
        dWin(i) = 1;
    elseif(discrete_result == 0 && continuous_result == 1)
        cWin(i) = 1;
    elseif(discrete_result == 1 && discrete_result == 1)
        tie(i) = 1;
    end
    
end