% Define the objective function
fun = @(x) [10*(x(2) - x(1)^2); 
            (1 - x(1))];

% Initial guess
x0 = [-1; 1];

% Lower and upper bounds
lb = [-inf; 0];
ub = [inf; inf];

% Options
options = optimoptions('lsqnonlin', 'Display', 'iter');

% Call lsqnonlin
[x, resnorm, residual, exitflag, output] = lsqnonlin(fun, x0, lb, ub, options);

% Display results
disp('Optimal solution:');
disp(x);
disp('Residual norm:');
disp(resnorm);
disp('Exit flag:');
disp(exitflag);
