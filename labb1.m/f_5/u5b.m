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

% Startparametrar
k1 = 7000; k2 = 100000;
guess = [k1, k2]';
k_0 = guess;

ck = 0.71; 
cs = 1;

% För att lagra värden och fel
k1_vals = [];
k2_vals = [];
errors = [];
accuracy_vals = [];

% Generera referenslösning (exakt lösning för jämförelse)
ref_k1 = 691.5829; % Exakt k1
ref_k2 = 182116.4403; % Exakt k2
ref_solution = [ref_k1; ref_k2];

% Newton-metoden för system
tolerance = 1e-6;
max_iteration = 100; 
accuracy = 10;
i = 0;

while accuracy > tolerance && i < max_iteration
    % Beräkna Jacobian och transferfunktioner
    J = Jacobian_transfer_functions(k1, k2);
    F = transfer_functions(k1, k2, ck, cs);
    b = J \ F;
    k_ny = k_0 - b; 
    accuracy = norm(k_ny - k_0);
    
    % Uppdatera variabler
    k_0 = k_ny;
    k1 = k_ny(1);  
    k2 = k_ny(2);  
    
    % Spara värden
    k1_vals = [k1_vals, k1];
    k2_vals = [k2_vals, k2];
    accuracy_vals = [accuracy_vals, accuracy];
    
    % Spara fel för konvergensstudie
    error = norm([k1 - ref_solution(1), k2 - ref_solution(2)]);
    errors = [errors, error];
    
    % Skriv ut resultat
    disp(['Iteration ', num2str(i), ' - k1: ', num2str(k1), ', k2: ', num2str(k2), ', Accuracy: ', num2str(accuracy)]);
    i = i + 1;
end

if accuracy <= tolerance
    disp('Konvergerat!');
else
    disp('Konvergerar inte.');
end

% Beräkning av konvergensordningen p
p_vals = [];
if length(errors) >= 3
    for j = 2:length(errors)-1
        p = log(errors(j+1)/errors(j)) / log(errors(j)/errors(j-1));
        disp(['p är ', num2str(p)]);
        p_vals = [p_vals, p];
    end
end

% Plottning av konvergens
figure;

% First subplot
subplot(3, 1, 1); % 3 rows, 1 column, first plot
plot(1:i, accuracy_vals, 'o-', 'DisplayName', 'Accuracy');
xlabel('Iteration');
ylabel('Accuracy');
title('Konvergensmätning');
grid on;
legend;

% Second subplot
subplot(3, 1, 2); % 3 rows, 1 column, second plot
plot(1:length(errors), errors, 'x-', 'DisplayName', 'Errors');
xlabel('Iteration');
ylabel('Error');
title('Fel');
grid on;
legend;

% Third subplot
subplot(3, 1, 3); % 3 rows, 1 column, third plot
if ~isempty(p_vals)
    plot(2:length(p_vals)+1, p_vals, 'o-', 'DisplayName', 'p (konvergensordning)');
    xlabel('Iteration');
    ylabel('Konvergensordning');
    title('Konvergensordning');
    grid on;
    legend;
end
