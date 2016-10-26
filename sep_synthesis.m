% simple synthesis test for 2 subsystems

A_ii = eye(2);
B_i = [1;1];

ss1 = Subsys(2, .01, [100 100], [0 0], .01, 130, -.65);
ss1.setd([0 0]);
ss1.setAB(A_ii, B_i);

ss1.getTrans();