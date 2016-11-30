% assumption mining for 2 subsystems

% subsystem dynamics
A_ii = [.1 .2; .3 .1];
A_ij = [.1 .1; .05 .05]; % = [.1 .05]' * [1 1]
B_i = [0; 1];
C_i = [1 1];
E_i = [.1; .05];

% subsystem 1
ss1 = Subsys(2, .1, [11 11], [0 0], .1, 21, -.65);
ss1.setAB(A_ii, B_i);

% subsystem 2
ss2 = Subsys(2, .1, [11 11], [0 0], .1, 21, -.65);
ss2.setAB(A_ii, B_i);

% search helper
sh = Search_Helper(.1, [10 10], [0 0]);
sh.set_If(E_i, E_i);

% search process
sample = sh.find_sample();

while(sample ~= -1)

    % get assumptions
    [d1, d2] = sh.get_assumptions(sample);
    
    % set assumptions for each system
    ss1.setd(d2);
    ss2.setd(d1);
    
    % require that each system meets assumptions
    ss1.setSafe(C_i, d1);
    ss2.setSafe(C_i, d2);
    ss1.updateInv();
    ss2.updateInv();
    
    % shrink down to invariant sets
    volume1 = ss1.ConInvOI();
    volume2 = ss2.ConInvOI();
    
    % result
    if(volume1 ~= 0 && volume2 ~= 0)
        result = 'success';
    else
        result = 'fail';
    end
    
    % use monotone specification ideas to update search
    sh.update_search(sample, result);
    
end