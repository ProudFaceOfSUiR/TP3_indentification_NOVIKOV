function J = jacobien(theta, p)

theta = theta(:);
p = p(:);
l_0 = p(1);
l1 = p(2);
l2 = p(3);
l3 = p(4);
t_1   = theta(1) + p(5);
t_2   = theta(2) + p(6);
t_3   = theta(3) + p(7);
xc = p(8); yc = p(9); ac = p(10);

t_12  = t_1 + t_2;
t_123 = t_12 + t_3;
c_1 = cos(t_1); s1 = sin(t_1);
c_12 = cos(t_12); s12 = sin(t_12);
c_123 = cos(t_123); s123 = sin(t_123);

x_0 = l1*c_1 + l2*c_12 + l3*c_123;
y_0 = p(1) + l1*s1 + l2*s12 + l3*s123;
dx = x_0 - xc;
dy = y_0 - yc;

ca = cos(ac); sa = sin(ac);
R = [ca sa; -sa ca];
xm = ca*dx + sa*dy;
ym = -sa*dx + ca*dy;

dX_dl0 = [0; 1];
dX_dl1 = [c_1; s1];
dX_dl2 = [c_12; s12];
dX_dl3 = [c_123; s123];

dX_do1 = [-(l1*s1 + l2*s12 + l3*s123);
          l1*c_1 + l2*c_12 + l3*c_123];
dX_do2 = [-(l2*s12 + l3*s123);
          l2*c_12 + l3*c_123];
dX_do3 = [-l3*s123;
          l3*c_123];

J = zeros(3,10);
J(1:2,1) = R*dX_dl0;
J(1:2,2) = R*dX_dl1;
J(1:2,3) = R*dX_dl2;
J(1:2,4) = R*dX_dl3;
J(1:2,5) = R*dX_do1;
J(1:2,6) = R*dX_do2;
J(1:2,7) = R*dX_do3;

J(1:2,8) = -R(:,1);
J(1:2,9) = -R(:,2);
J(1:2,10) = [ym; -xm];

J(3,5:7) = 1;
J(3,10)  = -1;
end
