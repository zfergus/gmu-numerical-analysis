% Computes the roots of a vector valued function using the Newton's Method.
% Written by Zachary Ferguson

function xc = multivariate_newtons_method(f, df, x0, tol, figHandle)
    % Compute the root to f(x) using Newton's Method
    % Input:
    %   f - vector valued function to find the roots of
    %   df - Jacobian of f(x)
    %   x0 - intial guess
    %   tol - tolerance for the root
    % Output:
    %   xc - computed root to the function f(x).
    if nargin < 4
        tol = 1e-8;
    end
    if nargin < 5
        figHandle = false;
    end

    n_steps = 0;
    xc = x0;
    fe = norm(f(xc), inf);
    errors = [fe];
    while fe > tol
        s = df(xc) \ -f(xc);
        xc = xc + s;
        n_steps = n_steps + 1;
        fe = norm(f(xc), inf);
        errors = [errors fe];
    end
    fprintf('\tNumber of steps to solve to %g accuracy: %d\n', tol, n_steps);

    if figHandle ~= false
        figure(figHandle);
        plot(1:size(errors, 2), errors, '-ob');
    end
end
