% Computes the fixed point of a function using the FPI.
% Written by Zachary Ferguson

function xc = fixed_point_iteration(g, x0, f, tol)
    % Compute the fixed point of g(x).
    % Input:
    %   g   - function to solve for the fixed point.
    %   x0  - initial guess
    %   f   - f(x) = g(x) - x
    %   tol - solution tolerance
    % Output:
    %   xc - computed root to the function g(x) = x.
    if nargin < 4
        tol = 1e-9;
    end

    r = fzero(f, x0);
    fprintf('r = %f\n', r);
    ei = 0;

    prev_x = x0;
    x = g(x0);
    n = 1;
    while (abs(prev_x - x) > 0.5 * tol)
        prev_x = x;
        x = g(x);
        n = n + 1;
        ei1 = abs(x - r);
        if (abs(prev_x - x) <= 0.5 * tol)
             fprintf('e_(i+1)/e_i = %.10f\n', ei1 / ei);
        end
        ei = ei1;
    end
    fprintf('n = %d\n', n);
    xc = x;
end
