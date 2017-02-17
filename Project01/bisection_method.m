% Computes the root of a function using the bisection method.
% Written by Zachary Ferguson

function r = bisection_method(f, a, b, tol)
    % Finds the root of the function, f, in the interval [a, b] within an
    % absolute tolerance.
    if nargin < 4
        tol = 1e-7;
    end
    assert(f(a) * f(b) <= 0);

    n = 0;
    while (abs(b - a) / 2.0) > (0.5 * tol)
        c = (a + b) / 2.0;
        if (f(c) == 0)
            break;
        end

        if (f(a) * f(c) <= 0)
            b = c;
        else
            a = c;
        end
        n = n + 1;
    end
    fprintf('\tn = %d\n', n);
    r = (a + b) / 2.0;
end
