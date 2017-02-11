% Computes the roots of a function using the Newton's Method.
% Written by Zachary Ferguson


% Solve the project questions
function newtons_method
    % Q1a
    fprintf('Q1a:\n\tf(x) = x^3 - 2x - 2 = 0');
    f  = @(x) x^3 - 2*x - 2;
    fp = @(x) 3*x^2 - 2;
    x0 = 2;
    fprintf('\tx0 = %.10f\n', x0);
    fprintf('\tr = %.10f\n', compute_root(f, fp, x0));

    % Q1b
    fprintf('Q1b:\n\tf(x) = e^x + x - 7 = 0');
    f  = @(x) exp(x) + x - 7;
    fp = @(x) exp(x) + 1;
    x0 = 0;
    fprintf('\tx0 = %.10f\n', x0);
    fprintf('\tr = %.10f\n', compute_root(f, fp, x0));

    % Q1c
    fprintf('Q1c:\n\tf(x) = e^x + sin(x) - 4 = 0');
    f  = @(x) exp(x) - sin(x) - 4;
    fp = @(x) exp(x) - cos(x);
    x0 = 2;
    fprintf('\tx0 = %.10f\n', x0);
    fprintf('\tr = %.10f\n', compute_root(f, fp, x0));
end


% Compute the foot to f(x)
function r = compute_root(f, fp, x0, tol)
    if nargin < 4
        tol = 1e-9;
    end

    n = 0;
    x = x0;
    while (abs(f(x)) > 0.5 * tol)
        x = x - f(x) / fp(x);
        n = n + 1
    end
    fprintf('\tn = %d\n', n)
    r = x
end
df(x0)
