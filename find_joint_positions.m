function [positions] = find_joint_positions(angles)

li = 15;
l1 = 50;
l2 = 41.5;
l3 = 55.8;
l4 = 40.1;
l5 = 39.4;
l6 = 61.9;
l7 = 39.3;
l8 = 36.7;
l9 = 49;

a = 7.8;
b = 38;

c = 1.729556;

ti = angles(1);
t1 = angles(2);
t2 = angles(3);
t3 = angles(4);
t4 = angles(5);
t5 = angles(6);
t6 = angles(7);
t7 = angles(8);
t8 = angles(9);

positions(1, 1)  = li*cos(ti);
positions(1, 2)  = li*sin(ti);
positions(2, 1)  = l1*cos(t1) + positions(1, 1);
positions(2, 2)  = l1*sin(t1) + positions(1, 2);
positions(3, 1)  = l3*cos(t3) + positions(2, 1);
positions(3, 2)  = l3*sin(t3) + positions(2, 2);
positions(4, 1)  = -b + l7*cos(t7);
positions(4, 2)  = -a + l7*sin(t7);
positions(5, 1)  = l8*cos(t8) + positions(4, 1);
positions(5, 2)  = l8*sin(t8) + positions(4, 2);
positions(6, 1)  = l9*cos(t8 + c) + positions(4, 1);
positions(6, 2)  = l9*sin(t8 + c) + positions(4, 2);
positions(7, 1) = 0;
positions(7, 2) = 0;
positions(8, 1) = -b;
positions(8, 2) = -a;

end
