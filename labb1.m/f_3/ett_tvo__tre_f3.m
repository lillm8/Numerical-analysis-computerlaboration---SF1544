m1=475;
m2=53;
k1 = 5400;
k2 = 135000;

c1=310;
c2=1200;
v=65/3.6;
H=0.24;
L=1;

% Starttid
t_0 = 0;

h = @(t) (t > L/v) .* 0 + (t <= L/v) .* ((H / 2) * (1 - cos((2 * pi * v * t) / L)));
dh_dt = @(t) (t > L/v) .* 0 + (t <= L/v) .* (H * pi * v / L) * sin((2 * pi * v * t) / L);

A = [0, 0, 1, 0;
  0, 0, 0, 1;
  -k1/m1, k1/m1, -c1/m1, c1/m1;
   k1/m2, -(k1 + k2)/m2, c1/m2, -(c1 + c2)/m2
];

g = @(t) [0; 0; 0; (c2 * dh_dt(t) + k2 * h(t)) / m2];


eigenvalues = eig(A);

%Ska bli vår lista med maximal tid%
ls = [];
for d = 1:length(eigenvalues)
    lambda = eigenvalues(d);
    ls(end + 1) = -2 * real(lambda) / abs(lambda)^2; 
end

delta_max = min(ls);

% Display result with formatted output
fprintf('Maximal tid blir %f\n', delta_max);

alpha = input('Korrigerande faktor alpha: ');

delta_max_ny = delta_max * alpha;

n=input('Ange antal steg: '); 
h_s=(delta_max_ny-t_0)/n; 
disp(['h= ' num2str(h_s) '.']); 


tspan=[t_0 delta_max_ny];
[tv,yv]=EulerSyst(@(t,y) f2by2(t, y, A, g(t)), tspan, y0, n);
plot(tv,yv(1,:),'b-', tv,yv(2,:), 'r-'); 
title('Eulersmetod tidsintervall %f, från alpa-faktor %f',delta_max_ny ,alpha);

