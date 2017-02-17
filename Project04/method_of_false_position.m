% Computes the roots of a function using the Method of False Position.
% Written by Zachary Ferguson

function xc = method_of_false_position(f, x0, x1, tol)
    % Compute the root to f(x) using the Method of False Position
    % Input:
    %   f - function to find the roots of
    %   x0, x1 - intial guess range (f(x0)*f(x1) < 0)
    %   tol - tolerance for the root
    % Output:
    %   xc - computed root to the function f(x).
    if nargin < 4
        tol = 1e-9;
    end

    n = 0;
    a = x0;
    b = x1;
    c = b;
    while (abs(f(c)) >= 0.5 * tol)
        if (f(a) * f(b)) < 0
            b = c;
        else
            a = c;
        end
        % Next value of c closer to the root
        c = (b * f(a) - a * f(b)) / (f(a) - f(b));
        n = n + 1;
    end
    fprintf('\tn = %d\n', n);
    xc = c;
end
