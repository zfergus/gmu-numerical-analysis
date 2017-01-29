% Computes the fixed point of a function using the FPI.
% Written by Zachary Ferguson


function fixed_point_iteration
    fprintf('Fixed Point Iteration\nWrtten by Zachary Ferguson\n');

    fprintf('Q1a:\n\tg(x) = (2x+2)^(1/3) = x');
    fprintf('\txc = %.10f', compute_fixed_point(lambda x: (2 * x + 2)**(1 / 3.),
        2))

    from math import log
    fprintf('Q1b:\n\tg(x) = ln(7-x) = x')
    fprintf('\txc = %.10f', compute_fixed_point(lambda x: log(7 - x), 2))

    from math import sin
    fprintf('Q1c:\n\tg(x) = ln(4-sin(x)) = x')
    fprintf('\txc = %.10f', compute_fixed_point(lambda x: log(4 - sin(x)), 2))

    A = 3.
    g = lambda x: (x + A / x) / 2.
    fprintf('Q3a:\n\tg(x) = (x + 3 / x) / 2')
    fprintf('\tx0 = 2')
    fprintf('\txc = %.10f', compute_fixed_point(g, 2, print_n=True))

    A = 5.
    fprintf('Q3b:\n\tg(x) = (x + 5 / x) / 2')
    fprintf('\tx0 = 2')
    fprintf('\txc = %.10f', compute_fixed_point(g, 2, print_n=True))

    from math import cos
    fprintf('Q5:\n\tg(x) = cos^2(x)')
    xc = compute_fixed_point(lambda x: (cos(x))**2, 1, tol=1e-6, print_n=True)
    fprintf('\txc = %.10f', xc)
    fprintf('\n\td/dx g(x) = -2*cos(x)*sin(x)')
    fprintf('\t|d\dx g(xc)| = %.10f', abs(-2 * cos(xc) * sin(xc)))
    fprintf('\tTherefore g(x) is locally convergent to xc.')
end


% Compute the fixed point of g(x).
function xc = compute_fixed_point(g, x0, tol)
    if nargin < 3
        tol = 1e-7;
    end

    prev_x = x0;
    x = g(x0);
    n = 1;
    while (abs(prev_x - x) > 0.5 * tol)
        prev_x = x;
        x = g(x);
        n += 1;
    end
    fprintf('\tn = %d', n);
    end
    xc = x;
end
