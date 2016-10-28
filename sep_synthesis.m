% simple synthesis test for 2 subsystems

A_ii = eye(2);
B_i = [1;1];

ss1 = Subsys(2, .2, [5 5], [0 0], .2, 10, -1);
ss1.setd([1 0; 0 1]);
ss1.setAB(A_ii, B_i);

ss1.getTrans();