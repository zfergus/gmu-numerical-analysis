% Computes the roots of a function using the Inverse Quadratic Interpolation.
% Written by Zachary Ferguson

function xc = inverse_quadratic_interpolation(f, x0, x1, x2, tol)
    % Compute the root to f(x) using the Inverse Quadratic Interpolation
    % Input:
    %   f - function to find the roots of
    %   x0, x1, x2 - intial guess defining the parabola
    %   tol - tolerance for the root
    % Output:
    %   xc - computed root to the function f(x).
    if nargin < 5
        tol = 1e-9;
    end

    n = 0;
    a = x0;
    b = x1;
    c = x2;
    div = @(x,y) (f(x) / f(y));
    while (abs(f(c)) >= 0.5 * tol)
        q = div(a, b);
        r = div(c, b);
        s = div(c, a);

        % Next value of c closer to the root
        d = c - (r * (r - q) * (c - b) + (1 - r) * s * (c - a)) / ...
            ((q - 1) * (r - 1) * (s - 1));

        a = b;
        b = c;
        c = d;
        n = n + 1;
    end
    fprintf('\tn = %d\n', n)
    xc = c;
end
