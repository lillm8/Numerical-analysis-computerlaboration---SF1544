m1=475;
m2=53;
k1 = 5400;
k2 = 135000;

c1=310;
c2=1200;
v=65/3.6;
H=0.24;
L=1;

% _____Euler frammåt_____


% Start och sluttid
t_0 = 0;
T_1 = 5*10^-4;
T_2 = 5*10^-3;
tspan_1=[t_0 T_1];
tspan_2=[t_0, T_2];
y0=[0; 0; 0; 0];


% Skapa funktionen h(t) och dess derivata
h = @(t) (t > L/v) .* 0 + (t <= L/v) .* ((H / 2) * (1 - cos((2 * pi * v * t) / L)));
dh_dt = @(t) (t > L/v) .* 0 + (t <= L/v) .* (H * pi * v / L) * sin((2 * pi * v * t) / L);


A = [0, 0, 1, 0;
  0, 0, 0, 1;
  -k1/m1, k1/m1, -c1/m1, c1/m1;
   k1/m2, -(k1 + k2)/m2, c1/m2, -(c1 + c2)/m2
];

g = @(t) [0; 0; 0; (c2 * dh_dt(t) + k2 * h(t)) / m2];

n=input('Ange antal steg: '); 
h_s=(T-t_0)/n; 
disp(['h= ' num2str(h_s) '.']); 

[tv,yv]=EulerSyst(@(t,y) f2by2(t, y, A, g(t)), tspan_1, y0, n);
plot(tv,yv(1,:),'b-', tv,yv(2,:), 'r-'); 
title('Eulersmetod tidsintervall %f, från alpa-faktor %f', , );