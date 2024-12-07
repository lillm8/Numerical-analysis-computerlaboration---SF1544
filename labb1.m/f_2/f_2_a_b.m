m1=475;
m2=53;
k1 = 5400;
k2 = 135000;

c1=310;
c2=1200;
v=65/3.6;
H=0.24;
L=1;

% Start och sluttid
t_0 = 0;
T = 1;

% Skapa funktionen h(t) och dess derivata
h = @(t) (t > L/v) .* 0 + (t <= L/v) .* ((H / 2) * (1 - cos((2 * pi * v * t) / L)));
dh_dt = @(t) (t > L/v) .* 0 + (t <= L/v) .* (H * pi * v / L) * sin((2 * pi * v * t) / L);

% Definiera ODE-systemet
% y(1) = z_1
% y(2) = z_2
% y(3) = dotderiv z_1
% y(4) = dotderiv z_2

A = [0, 0, 1, 0;
  0, 0, 0, 1;
  -k1/m1, k1/m1, -c1/m1, c1/m1;
   k1/m2, -(k1 + k2)/m2, c1/m2, -(c1 + c2)/m2
];

g = @(t) [0; 0; 0; (c2 * dh_dt(t) + k2 * h(t)) / m2];


ode_system = @(t,y) A * y + g(t);


% Initialvärden
y0 = [0; 0; 0; 0]; % [y1_dot, y2_dot, y_1, y_2] NOTERA årdningen av variabeln ode_system.

t_span = [t_0, T];    % tidsintervall

% RelTol sätter det maximala tolererade relativa felet i varje steg. 
% Det finns ingen garanti att det relativa globala felet i en viss tidpunkt ligger under detta värde.
% AbsTol talar om det värde på storleken på lösningen under vilket det
% relativa felet inte ska kontrolleras.

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

% Lös ODE-systemet MED implementation av "options".
[t, y] = ode45(ode_system, t_span, y0, options);

% Plotta resultaten
plot(t, y(:, 1), '-b', t, y(:, 2), '-r');
legend('z_1','z_2')

max_z1 = max(y(: , 1));
max_z2 = max(y(: , 2));

fprintf('Maximalt utslag för z_1 är %f och maximalt utslag för z_2 är %f', max_z1, max_z2);

