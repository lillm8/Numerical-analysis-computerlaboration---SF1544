addpath('./');

% Parameterinställningar
m1 = 475;
m2 = 53;
k1 = 5400;
k2 = 135000;
c1 = 310;
c2 = 1200;
v = 65 / 3.6;
H = 0.24;
L = 1;

% Initialvärden
t_0 = 0;
T = 3;
y0 = [0; 0; 0; 0]; % [z1, z2, z1_dot, z2_dot]

% Funktioner
h = @(t) (t > L/v) .* 0 + (t <= L/v) .* ((H / 2) * (1 - cos((2 * pi * v * t) / L)));
dh_dt = @(t) (t > L/v) .* 0 + (t <= L/v) .* (H * pi * v / L) * sin((2 * pi * v * t) / L);

A = [0, 0, 1, 0;
     0, 0, 0, 1;
    -k1/m1, k1/m1, -c1/m1, c1/m1;
     k1/m2, -(k1 + k2)/m2, c1/m2, -(c1 + c2)/m2];

g = @(t) [0; 0; 0; (c2 * dh_dt(t) + k2 * h(t)) / m2];

% Egenvärden och delta_max
ev = eig(A);
ls = arrayfun(@(lambda_k) -2 / (real(lambda_k)^2 + imag(lambda_k)^2) * real(lambda_k), ev);
delta_max = min(ls);
fprintf('Maximal tid: %f\n', delta_max);

% Värden av alpha och antal steg
alpha_values = [0.9, 1.0, 1.1, 1.5];

% Initialisera figur
figure('Position', [100, 100, 1200, 800]); % Större fönsterstorlek

% Kör Euler-metoden för varje alpha
for i = 1:length(alpha_values)
    alpha = alpha_values(i);
    delta_max_ny = delta_max * alpha;
    h_s = delta_max_ny; % Steglängd
    n = round((T-t_0)/h_s);
    fprintf('Alpha = %f, Δt = %f, h = %f\n', alpha, delta_max_ny, h_s);
    
    tspan = [t_0, T];
    [tv, yv] = EulerSyst(@(t, y) f2by2(y, A, g(t)), tspan, y0, n);
    
    % Skapa subplot för denna alpha
    subplot(2, 2, i); % 2x2 layout, aktuell subplot är i
    plot(tv, yv(1, :), '-b', 'DisplayName', 'z_1');
    hold on;
    plot(tv, yv(2, :), '-r', 'DisplayName', 'z_2');
    hold off;
    
    % Anpassa subplot
    title(sprintf('Euler-metoden: α = %.1f', alpha));
    xlabel('Tid (s)');
    ylabel('Position');
    legend('show', 'Location', 'best');
    grid on;
end

% Justera layout
sgtitle('Euler-metoden för olika värden av \alpha');
