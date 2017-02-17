% Computes the roots of a function using the Secant Method.
% Written by Zachary Ferguson

function xc = secant_method(f, x0, x1, tol)
    % Compute the root to f(x) using the Secant Method
    % Input:
    %   f - function to find the roots of
    %   x0, x1 - intial guesses
    %   tol - tolerance for the root
    % Output:
    %   xc - computed root to the function f(x).
    if nargin < 4
        tol = 1e-9;
    end

    n = 0;
    xi = x1;
    xi_1 = x0;
    while (abs(f(xi)) >= 0.5 * tol)
        x = xi - (f(xi) * (xi - xi_1)) / (f(xi) - f(xi_1));
        xi_1 = xi;
        xi = x;
        n = n + 1;
    end
    fprintf('\tn = %d\n', n)
    xc = xi;
end
