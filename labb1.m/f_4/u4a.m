% Kod som plotta alpha = 1, = 10, = 100 bredvid varandra
% Det som händer stämmer med teorin för den implicita trapetzmetoden
% Det vi ser är att alpha = 1 och alpha = 10 är väldigt lika varandra trots
% den stora skillnaden i steghöjden, vilket är förväntat.

addpath('./');
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

% Initial conditions
t_0 = 0;
T = 3;
y0 = [0; 0; 0; 0];

h = @(t) (t > L/v) .* 0 + (t <= L/v) .* ((H / 2) * (1 - cos((2 * pi * v * t) / L)));
dh_dt = @(t) (t > L/v) .* 0 + (t <= L/v) .* (H * pi * v / L) * sin((2 * pi * v * t) / L);

A = [0, 0, 1, 0;
     0, 0, 0, 1;
    -k1/m1, k1/m1, -c1/m1, c1/m1;
     k1/m2, -(k1 + k2)/m2, c1/m2, -(c1 + c2)/m2
];

g = @(t) [0; 0; 0; (c2 * dh_dt(t) + k2 * h(t)) / m2];

ev = eig(A);

ls = [];

for k = 1:length(ev)
    lambda_k = ev(k);
    Re_lambda_k = real(lambda_k);
    abs_lambda_k_sq = real(lambda_k)^2 + imag(lambda_k)^2;
    ls(k) = -2 / abs_lambda_k_sq * Re_lambda_k;
end

delta_max = min(ls);
fprintf('The maximal time is %f\n', delta_max);

alphas = [1, 10, 100];

figure;
for i = 1:length(alphas)
    alpha = alphas(i);
    delta_max_ny = delta_max * alpha;
    h_s = delta_max_ny;
    n = round((T - t_0) / h_s);
    
    disp(['h= ' num2str(h_s) '.']); 
    theta = 0.5;

    tspan = [t_0, T];
    [tv, yv] = EulerSyst_4(@(t, y) f2by2(y, A, g(t)), tspan, y0, n, theta, A, h_s, g);
    
    % Plot in the appropriate subplot
    subplot(1, 3, i);
    plot(tv, yv(1, :), 'b-', tv, yv(2, :), 'r-');
    title(sprintf('α = %d', alpha));
    xlabel('Time (s)');
    ylabel('Displacement');
    legend('c(1)', 'c(2)');
end
