% Initial guess

% Tre gissningar som konvergerar till rätt k_1:
% (k1, k2) = 
% (3000, 240000)
% (2000, 80000)
% (7000, 100000)


% Tre gissningar som inte konvergerar: 
% (k1, k2) = 
% (70000, 5000)
% (400000, 80000)
% (700000, 15400)



k1 = 3000; k2 = 24000;
guess = [k1, k2]';
k_0 = guess;

ck = 0.71; 
cs = 1;

% För att lagra k1 och k2 under iterationer
k1_vals = [];
k2_vals = [];

% The Newton method for systems 
accuracy = 10;
tolerance = 1e-6;
i = 0;
max_iteration = 100; 

while accuracy > tolerance && i < max_iteration
    J = Jacobian_transfer_functions(k1, k2);
    F = transfer_functions(k1, k2, ck, cs);
    b = J \ F;
    k_ny = k_0 - b; 
    accuracy = norm(k_ny - k_0);
    k_0 = k_ny;
    k1 = k_ny(1);  
    k2 = k_ny(2);  
    
    % Spara värden
    k1_vals = [k1_vals, k1];
    k2_vals = [k2_vals, k2];
    
    disp([' k1: ', num2str(k1), ' k2: ', num2str(k2)]);
    i = i + 1;
end

if accuracy <= tolerance
    disp('Kovergerat!');
else
    disp('Konververat inte.');
end

% Plottning
figure;
iter = 1:length(k1_vals); % Iterationsnummer
plot(k2_vals, iter, 'x-', 'DisplayName', 'k2');

xlabel('Values of  k2');
ylabel('Iteration');
title('Convergence of k2');
legend;
grid on;
hold off;

