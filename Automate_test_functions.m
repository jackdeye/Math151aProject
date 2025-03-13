% Define functions and analytical solutions
functions = {@(x) 1./(1+1000*(x-0.5).^2), ...      % f1
             @(x) 1./x.*sqrt((1+x)./(1-x)).*log((2*x.^2+2*x+1)./(2*x.^2-2*x+1)), ... % f2
             @(x) sqrt(abs(x)), ...                % f3
             @(x) 1./(1+x.^2), ...                 % f4
             @(x) cos(x)};                         % f5
solutions = [atan(5*sqrt(10))/(5*sqrt(10)), ...    % e1
             4*pi*acot(sqrt((1 + sqrt(5))/2)), ... % e2
             4/3, ...                              % e3
             pi/4, ...                             % e4
             2];                                   % e5
bounds = {[0, 1], ...                              % Bounds for f1
          [-1, 1], ...                             % Bounds for f2
          [-1, 1], ...                             % Bounds for f3
          [0, 1], ...                              % Bounds for f4
          [-pi/2, pi/2]};                          % Bounds for f5

% Create output folder
output_folder = 'data';
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Number of iterations and skipping steps
num = 80;
skip = 2;
begin = 2;
x = 1:skip:num;

methods = {'GQ', 'CS', 'CT'}; % Gauss-Legendre, Composite Simpson's, Composite Trapezoid

% Loop through all functions
for i = 1:length(functions)
    f = functions{i};          
    exact_solution = solutions(i); 
    bound = bounds{i};

    % Define evaluation functions for each method
    eval_methods = {
        @(n) gauss_legendre_quadrature(f, bound(1), bound(2), n), ...
        @(n) composite_simpsons_rule(f, bound(1), bound(2), n), ...
        @(n) composite_trapezoid_rule(f, bound(1), bound(2), n)
    };

    for j = 1:length(eval_methods)
        eval = eval_methods{j}; 
        values = zeros(num/skip, 1); 
        for n = begin:skip:num
            diff = exact_solution - eval(n);
            values(n/skip) = diff;
            disp(['Function f' num2str(i) ', Method ' methods{j} ', n = ' num2str(n)]);
            if abs(diff) < eps
                break;
            end
        end

        % Plot and save linear scale graph
        figure;
        plot(x, values(:), 'Color', 'magenta');
        title(['Function f' num2str(i) ' - Linear Scale (' methods{j} ')']);
        xlabel('n');
        ylabel('Error');
        saveas(gcf, fullfile(output_folder, [methods{j} '-f' num2str(i) '-linspace.fig']));

        % Plot and save log-log scale graph
        figure;
        plot(log(x), log(values(:)), 'Color', 'green');
        title(['Function f' num2str(i) ' - Log-Log Scale (' methods{j} ')']);
        xlabel('log(n)');
        ylabel('log(Error)');
        saveas(gcf, fullfile(output_folder, [methods{j} '-f' num2str(i) '-logspace.fig']));
    end
end
