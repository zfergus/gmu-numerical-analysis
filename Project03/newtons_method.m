% Computes the roots of a function using the Newton's Method.
% Written by Zachary Ferguson

function xc = newtons_method(f, fp, x0, tol, m, print_ei)
    % Compute the root to f(x) using Newton's Method
    % Input:
    %   f - function to find the roots of
    %   fp - first dirivative of f(x)
    %   x0 - intial guess
    %   tol - tolerance for the root
    %   m - multiplicity of the root
    %   print_ei - which e_i limit should be printed
    % Output:
    %   xc - computed root to the function f(x).
    if nargin < 4
        tol = 1e-9;
    end
    if nargin < 5
        m = 1;
    end
    if nargin < 6
        print_ei = 0;
    end

    r = fzero(f, x0);

    n = 0;
    x = x0;
    ei = 1;
    ei_1 = 1;
    while (abs(f(x)) >= 0.5 * tol)
        x = x - m * f(x) / fp(x);
        n = n + 1;
        ei_1 = ei;
        ei = abs(r - x);
        if print_ei == 1
            fprintf('\te_i = %.8f; e_(i+1)/e_i = %.8f\n', ei, ei/ei_1);
        elseif print_ei == 2
            fprintf('\te_i = %.8f; e_(i+1)/(e_i)^2 = %.8f\n', ei, ei/(ei_1^2));
        end
    end
    fprintf('\tn = %d\n', n)
    xc = x;
end
