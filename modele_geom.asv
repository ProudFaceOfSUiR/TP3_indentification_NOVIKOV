function X = modele_geom(theta, p)

theta = theta(:);
p = p(:);
l_0 = p(1);
l1 = p(2);
l2 = p(3);
l3 = p(4);
t_1 = theta(1) + p(5);
t_2 = theta(2) + p(6);
t_3 = theta(3) + p(7);
xc = p(8);
yc = p(9);
ac = p(10);

t_12 = t_1 + t_2;
t_123 = t_12 + t_3;

x0 = l1*cos(t_1) + l2*cos(t_12) + l3*cos(t_123);
y0 = l_0 + l1*sin(t_1) + l2*sin(t_12) + l3*sin(t_123);
a0 = t_123 - pi/2;

dx = x0 - xc;  dy = y0 - yc;
xm = cos(ac)*dx + sin(ac)*dy;
ym = -sin(ac)*dx + cos(ac)*dy;
am = a0 - ac;

X = [xm; ym; am];
end
