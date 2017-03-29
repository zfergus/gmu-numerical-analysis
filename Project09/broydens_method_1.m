% Computes the roots of a vector valued function using the Broden's Method I.
% Written by Zachary Ferguson

function xc = broydens_method_1(f, A0, x0, tol, figHandle)
    % Compute the root to f(x) using Newton's Method
    % Input:
    %   f - vector valued function to find the roots of
    %   A0 - inital approximation for the Jacobian of f(x)
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
    A = A0;

    fe = norm(f(xc), inf);
    errors = [fe];
    while fe > tol
        s = A \ -f(xc);
        x_prev = xc;
        xc = xc + s;
        delta_f = f(xc) - f(x_prev);
        delta_x = xc - x_prev;
        A = A + ((delta_f - A * delta_x) * delta_x') / (delta_x' * delta_x);

        n_steps = n_steps + 1;
        fe = norm(f(xc), inf);
        errors = [errors fe];
    end
    fprintf('\tNumber of steps to solve to %g accuracy: %d\n', tol, n_steps);

    % Display the errors per iteration.
    if figHandle ~= false
        figure(figHandle);
        plot(1:size(errors, 2), errors, '-xr');
    end
end
