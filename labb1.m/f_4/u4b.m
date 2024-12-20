addpath('./');

% Parametrar
m1 = 475;
m2 = 53;
k1 = 5400;
k2 = 135000;
k2 = 100 * k2;

c1 = 310;
c2 = 1200;
v = 65 / 3.6;
H = 0.24;
L = 1;

% Initialvärden
t_0 = 0;
T = 0.05;
y0 = [0; 0; 0; 0];

% Funktioner
h = @(t) (t > L/v) .* 0 + (t <= L/v) .* ((H / 2) * (1 - cos((2 * pi * v * t) / L)));
dh_dt = @(t) (t > L/v) .* 0 + (t <= L/v) .* (H * pi * v / L) * sin((2 * pi * v * t) / L);

A = [0, 0, 1, 0;
     0, 0, 0, 1;
     -k1/m1, k1/m1, -c1/m1, c1/m1;
      k1/m2, -(k1 + k2)/m2, c1/m2, -(c1 + c2)/m2];

g = @(t) [0; 0; 0; (c2 * dh_dt(t) + k2 * h(t)) / m2];

% Egenvärden och maximalt tidssteg
ev = eig(A);
ls = [];
for k = 1:length(ev)
    lambda_k = ev(k);
    Re_lambda_k = real(lambda_k);
    abs_lambda_k_sq = real(lambda_k)^2 + imag(lambda_k)^2;
    ls(k) = -2 / abs_lambda_k_sq * Re_lambda_k;
end
delta_max = 1.25e-4; % min(ls)
fprintf('Maximalt tidssteg: %f\n', delta_max);

% Generera referenslösning med ode45
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
[ref_t, ref_y] = ode45(@(t, y) f2by2(y, A, g(t)), [t_0 T], y0, options);

% Konvergensstudie
theta = 0.5; % Trapetsmetoden
delta_t_vals = delta_max * [1, 1/2, 1/4, 1/8];
max_errors = [];

ero = 1;
i = 0;
for dt = delta_t_vals
    h_s = dt;
    n = round((T - t_0) / h_s);
    [tv, yv] = EulerSyst_4(@(t, y) f2by2(y, A, g(t)), [t_0 T], y0, n, theta, A, h_s, g);
    
    % Interpolera till referenslösningens tidspunkter
    yv_interp = interp1(tv, yv', ref_t)';

    % Beräkna max abs-fel i komponent z_2
    error1 = max(abs(yv_interp(1, :) - ref_y(:, 1)'));
    error2 = max(abs(yv_interp(2, :) - ref_y(:, 2)'));

    
    fprintf('Max error c_1: %.15f\n', error1);
    fprintf('Max error c_2: %.15f\n', error2);
     
    if i == 1;
        fprintf('Konvergensordningen är %f\n', log(ero/error2)/log(2));
    end
    i = 1;
    ero = error2;
   
end

plot(ref_t, ref_y(:, 2), 'k-', 'DisplayName', 'Referenslösning');
hold on;
plot(ref_t, yv_interp(2, :), 'r--', 'DisplayName', 'EulerSyst\_4');
legend;
hold off;





