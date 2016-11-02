% simple synthesis test for 2 subsystems

A_ii = eye(2);
B_i = [1;1];

ss1 = Subsys(2, .8, [6 6], [0 0], .2, 11, -1);
ss1.setd([0 0; 1 1]);
ss1.setAB(A_ii, B_i);

ss1.getTrans();

ss1.setGen(2);
ss1.setInv([2 4; 4 2]);
result = ss1.verifyInv();