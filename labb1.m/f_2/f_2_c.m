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
T_1 = 5*10^-4;
T_2 = 5*10^-3;
tspan_1=[t_0, T_1];
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



%_____ODE-45 delen_____
ode_system = @(t,y) A * y + g(t);
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);


% Sammanställt output av ODE-45 och Eulers metod
[t, y] = ode45(ode_system, tspan_2, y0, options);
[tv_1,yv_1]=EulerSyst(@(t,y) f2by2(y, A, g(t)), tspan_1, y0, n);
[tv_2,yv_2]=EulerSyst(@(t,y) f2by2(y, A, g(t)), tspan_2, y0, n);


figure('Position',[100, 100, 800, 800]);

% ODE-45
axes('Position', [0.1, 0.55, 0.8, 0.4]);
plot(t, y(:,1), 'r-', t, y(:,2), 'b-'); 
title('ODE-45 inom tidsintervallet 0 till 5*10^-3');
xlabel('Tid');
ylabel('Position');
legend('z_1','z_2')


% Euler med 10^-3 och 10^-4
axes('Position', [0.1, 0.06, 0.8, 0.4]); % [x, y, width, height]
plot(tv_1,yv_1(1,:),'b-', tv_1,yv_1(2,:),'b--', tv_2,yv_2(1,:),'r-', tv_2,yv_2(2,:),'r--'); 
title('Eulersmetod tidsintervall 5*10^-4 och 5*10^-3');
xlabel('Tid');
ylabel('Position');
legend('z_1 för 5*10^-3','z_2 för 5*10^-3', 'z_1 för 5*10^-4','z_1 för 5*10^-4')



