% simple synthesis test for 2 subsystems

% subsystem dynamics
A_ii = [.1 .2; .3 .1];
A_ij = [0 .3; .3 0];
B_i = [0; 1];

% subsystem 1
ss1 = Subsys(2, .1, [21 21], [-1 -1], .1, 21, -.65);
ss1.setAB(A_ii, B_i);
ss1.setGen(3);

% subsystem 2
ss2 = Subsys(2, .1, [21 21], [-1 -1], .1, 21, -.65);
ss2.setAB(A_ii, B_i);
ss2.setGen(3);

% choosing generator sample points
gen_offset1 = [5 10; 4 8; 3 6; 2 4; 1 2];
gen_offset2 = [10 10; 8 8; 6 6; 4 4; 2 2];
gen_offset3 = [10 5; 8 4; 6 3; 4 2; 2 1];

% add in offsets
gen_upper1 = 10*ones(size(gen_offset1)) + gen_offset1;
gen_lower1 = 10*ones(size(gen_offset1)) - gen_offset1;
gen_upper2 = 10*ones(size(gen_offset2)) + gen_offset2;
gen_lower2 = 10*ones(size(gen_offset2)) - gen_offset2;
gen_upper3 = 10*ones(size(gen_offset3)) + gen_offset3;
gen_lower3 = 10*ones(size(gen_offset3)) - gen_offset3;

% convert to state coordinates
x_upper1 = gen_upper1.*.1 - ones(size(gen_upper1)) + .05.*ones(size(gen_upper1));
x_lower1 = gen_lower1.*.1 - ones(size(gen_lower1)) - .05.*ones(size(gen_lower1));
x_upper2 = gen_upper2.*.1 - ones(size(gen_upper2)) + .05.*ones(size(gen_upper2));
x_lower2 = gen_lower2.*.1 - ones(size(gen_lower2)) - .05.*ones(size(gen_lower2));
x_upper3 = gen_upper3.*.1 - ones(size(gen_upper3)) + .05.*ones(size(gen_upper3));
x_lower3 = gen_lower3.*.1 - ones(size(gen_lower3)) - .05.*ones(size(gen_lower3));

% brute force search for invariant sets
check = 0;
best = 0;
inv1 = [];
inv2 = [];
for i = 1:5
    for j = 1:5
        for k = 1:5
            for l = 1:5
                for m = 1:5
                    for n = 1:5
                        
                        % move generators
                        upper1 = [gen_upper1(i, :); gen_upper2(j, :); gen_upper3(k, :)];
                        lower1 = [gen_lower1(i, :); gen_lower2(j, :); gen_lower3(k, :)];
                        upper2 = [gen_upper1(l, :); gen_upper2(m, :); gen_upper3(n, :)];
                        lower2 = [gen_lower1(l, :); gen_lower2(m, :); gen_lower3(n, :)];
                        
                        % if this is a subset, don't bother checking it
                        if(best ~= 0 && sum(sum(upper1 <= inv1)) == 6 && sum(sum(upper2 <= inv2)))
                            continue;
                        end
                        
                        % associated state coordinates
                        xup1 = [x_upper1(i, :); x_upper2(j, :); x_upper3(k, :)];
                        xlo1 = [x_lower1(i, :); x_lower2(j, :); x_lower3(k, :)];
                        xup2 = [x_upper1(l, :); x_upper2(m, :); x_upper3(n, :)];
                        xlo2 = [x_lower1(l, :); x_lower2(m, :); x_lower3(n, :)];
                        
                        % update invariant set, disturbances, transition
                        ss1.setInv(upper1, lower1);
                        ss1.setd([(A_ij*xlo1')'; (A_ij*xup1')']);
                        ss1.getTrans();
                        ss2.setInv(upper2, lower2);
                        ss2.setd([(A_ij*xlo2')'; (A_ij*xup2')']);
                        ss2.getTrans();
                        
                        % verify invariant set
                        result = ss1.verifyInv() + ss2.verifyInv();
                        if(result == 2)
                            vol = ss1.volume() + ss2.volume();
                            if(vol > best)
                                best = vol;
                                inv1 = upper1;
                                inv2 = upper2;
                            end
                        end
                        
                        check = check + 1
                    end
                end
            end
        end
    end
end