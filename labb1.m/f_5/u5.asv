% Initial guess
k1 = 5400; k2 = 135000;
guess  = [k1, k2]';
k_0 = guess;  

ck = 0.71; 
cs = 1;

% The Newton method for systems 
accuracy = 10;
tolerance = 1e-6;
i = 0;
max_iteration = 100; 

while accuracy > tolerance && i < max_iteration
    J = Jacobian_transfer_functions(k1, k2);
    F = transfer_functions (k1, k2, ck, cs);
    b = J \ F;
    k_ny = k_0 - b; 
    accuracy = norm(k_ny - k_0);
    k_0 = k_ny;
    k1 = k_ny(1);  
    k2 = k_ny(2);  
    disp([' k1: ', num2str(k1), ' k2: ', num2str(k2)]);
    i = i + 1;
end

if accuracy <= tolerance
    disp('Kovergerat!');
else
    disp('Konververat inte.');
end
