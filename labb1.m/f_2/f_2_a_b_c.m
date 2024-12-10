addpath('./');

m1=475;
m2=53;
k1 = 5400;
k2 = 135000;

c1=310;
c2=1200;
v=65/3.6;
H=0.24;
L=1;



% _____Eulerdelen_____


% Start och sluttid
t_0 = 0;
T = 1;
tspan=[t_0, T];
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


h1 = 5*10^(-3); 
h2 = 5*10^(-4);

n1 = round((T - t_0)/h1);
n2 = round((T - t_0)/h2);


%_____ODE-45 delen_____
ode_system = @(t,y) A * y + g(t);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);


% Sammanställt output av ODE-45 och Eulers metod
[t, y] = ode45(ode_system, tspan, y0, options);
[tv_1,yv_1]=EulerSyst(@(t,y) f2by2(y, A, g(t)), tspan, y0, n1);
[tv_2,yv_2]=EulerSyst(@(t,y) f2by2(y, A, g(t)), tspan, y0, n2);


figure('Position',[100, 100, 800, 800]);

% ODE-45
axes('Position', [0.1, 0.55, 0.8, 0.4]);
plot(t, y(:,1), 'r--', t, y(:,2), 'r',tv_1,yv_1(1,:),'g--', tv_1,yv_1(2,:),'g', tv_2,yv_2(1,:),'b--', tv_2,yv_2(2,:),'b' ); 
title('ODE-45 + Eulersmetod tidssteg 5*10^-4 och 5*10^-3');
xlabel('Tid');
ylabel('Position');
legend('z_1 ODE','z_2 ODE', 'z_1 EUL, 5*10^-3', 'z_2 EUL, 5*10^-3', 'z_1 EUL, 5*10^-4', 'z_1 EUL, 5*10^-4')


% Euler med 10^-3 och 10^-4.
axes('Position', [0.1, 0.06, 0.8, 0.4]); % [x, y, width, height]
plot(tv_1,yv_1(1,:),'g--', tv_1,yv_1(2,:),'g', tv_2,yv_2(1,:),'b--', tv_2,yv_2(2,:),'b');
title('Eulersmetod tidssteg 5*10^-4 och 5*10^-3');
xlabel('Tid');
ylabel('Position');
legend('z_1 EUL, 5*10^-3', 'z_2 EUL, 5*10^-3', 'z_1 EUL, 5*10^-4', 'z_1 EUL, 5*10^-4')


max_z1 = max(y(: , 1));
max_z2 = max(y(: , 2));

fprintf('Maximalt utslag för z_1 är %f och maximalt utslag för z_2 är %f', max_z1, max_z2);
