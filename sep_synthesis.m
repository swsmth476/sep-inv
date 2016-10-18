% simple synthesis test for 2 subsystems

% each block matrix has negative eigenvalues
A_ii = [.1 .2; .3 .1];
A_ij = [0 .3; .3 0];
B_i = [1; 0];
E_i = [1; 1];

ss1 = Subsys(2, .01, [100 100]', [0 0]', .01, 65, -.65);