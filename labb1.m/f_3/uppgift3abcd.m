addpath('./');
m1=475;
m2=53;
k1 = 5400;
k2 = 135000;

is_question_d = true;

% Question (d)
% k2 = 100 * k2
if is_question_d
    k2 = 100 * k2;
end

c1=310;
c2=1200;
v=65/3.6;
H=0.24;
L=1;

% Initial conditions
t_0 = 0;
T = 1;
y0 = [0; 0; 0; 0];

h = @(t) (t > L/v) .* 0 + (t <= L/v) .* ((H / 2) * (1 - cos((2 * pi * v * t) / L)));
dh_dt = @(t) (t > L/v) .* 0 + (t <= L/v) .* (H * pi * v / L) * sin((2 * pi * v * t) / L);

A = [0, 0, 1, 0;
  0, 0, 0, 1;
  -k1/m1, k1/m1, -c1/m1, c1/m1;
   k1/m2, -(k1 + k2)/m2, c1/m2, -(c1 + c2)/m2
];

g = @(t) [0; 0; 0; (c2 * dh_dt(t) + k2 * h(t)) / m2];

% Find all the eigenvalues of A
ev = eig(A);

% ls is the list containing the "maximal" time for each lambda
ls = [];

for k = 1:length(ev)
    lambda_k = ev(k);
    Re_lambda_k = real(lambda_k);
    abs_lambda_k_sq = real(lambda_k)^2 + imag(lambda_k)^2;
    ls(k) = -2 / abs_lambda_k_sq * Re_lambda_k;
end

% finds the minimal maximal eigenvalue, which gives the minimal time (b)
delta_max = min(ls);

% Display result with formatted output
fprintf('The maximal time is %f\n', delta_max);

if is_question_d
    alpha = 0.1;
else
    alpha = input('Corrected alpha factor: ');
end


delta_max_ny = delta_max * alpha;

% \Delta max
% d) - asks us to change alpha = 0.1 and observe the new \Delta_max, which
% is equal to 1.1181e-05
disp(delta_max_ny)

h_s=delta_max_ny; 
n = round((T-t_0)/h_s);
disp(['h= ' num2str(h_s) '.']);

tspan=[t_0, T];
[tv,yv]=EulerSyst(@(t,y) f2by2(y, A, g(t)), tspan, y0, n);
plot(tv,yv(1,:),'b-', tv,yv(2,:), 'r-'); 
title(sprintf('Euler method: Δt=%f, α=%f, n=%f, ', delta_max_ny, alpha, n));