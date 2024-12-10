addpath('./');

m1 = 475;
m2 = 53;

% Referensvärden
k1_ref = 5400;
k2_ref = 135000;

% Faktiska värden
k1 = 691.5829;
k2 = 182116.4403;

c1 = 310;
c2 = 1200;
v = 65 / 3.6;
H = 0.24;
L = 1;

% Start och sluttid
t_0 = 0;
T = 1;

% Skapa funktionen h(t) och dess derivata
h = @(t) (t > L/v) .* 0 + (t <= L/v) .* ((H / 2) * (1 - cos((2 * pi * v * t) / L)));
dh_dt = @(t) (t > L/v) .* 0 + (t <= L/v) .* (H * pi * v / L) * sin((2 * pi * v * t) / L);

% ODE-systemets matris A och funktion g(t)
ode_system = @(t, y, k1, k2) [
    0, 0, 1, 0;
    0, 0, 0, 1;
    -k1/m1, k1/m1, -c1/m1, c1/m1;
    k1/m2, -(k1 + k2)/m2, c1/m2, -(c1 + c2)/m2
] * y + [0; 0; 0; (c2 * dh_dt(t) + k2 * h(t)) / m2];

% Initialvärden
y0 = [0; 0; 0; 0]; % [y_1, y_2, dot(y_1), dot(y_2)]
t_span = [t_0, T]; % tidsintervall

% Lös ODE-systemet för aktuella värden
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t, y] = ode45(@(t, y) ode_system(t, y, k1, k2), t_span, y0, options);

% Lös ODE-systemet för referensvärden
[t_ref, y_ref] = ode45(@(t, y) ode_system(t, y, k1_ref, k2_ref), t_span, y0, options);

% Plotta resultaten
figure;
plot(t, y(:, 1), '-b', 'DisplayName', 'z_1 (actual)');
hold on;
plot(t, y(:, 2), '-r', 'DisplayName', 'z_2 (actual)');
plot(t_ref, y_ref(:, 1), '--b', 'DisplayName', 'z_1 (reference)');
plot(t_ref, y_ref(:, 2), '--r', 'DisplayName', 'z_2 (reference)');

title("Lösning av ODE med aktuella och referensvärden för k_1 och k_2");
xlabel('Tid (s)');
ylabel('Förskjutning (m)');
legend;
grid on;

% Maximala utslag
max_z1 = max(y(:, 1));
max_z2 = max(y(:, 2));
max_z1_ref = max(y_ref(:, 1));
max_z2_ref = max(y_ref(:, 2));

fprintf('Maximalt utslag för z_1 (actual) är %f och för z_2 (actual) är %f\n', max_z1, max_z2);
fprintf('Maximalt utslag för z_1 (reference) är %f och för z_2 (reference) är %f\n', max_z1_ref, max_z2_ref);