% Computes the fixed point of a function using the FPI.
% Written by Zachary Ferguson


function fixed_point_iteration
    fprintf('Fixed Point Iteration\nWrtten by Zachary Ferguson\n\n');

    fprintf('f(x) = 3*x^3 - 7*x^2 + 3*x - e^x + 2 = 0\n\n');
    f = @(x) 3*x^3 - 7*x^2 + 3*x - exp(x) + 2;
    
    fprintf('g1(x) = (e^x - 2) / (3x^2 - 7x + 3) = x\n')
    g1 = @(x) (exp(x) - 2) / (3*x^2 - 7*x + 3);
    fprintf('r1 = %.10f\n\n', compute_fixed_point(g1, -1, f));
    
    fprintf('g2(x) = (3*x^4 - 7*x^3 + 3*x^2 + 2*x) / e^x = x\n');
    g2 = @(x) (3*x^4 - 7*x^3 + 3*x^2 + 2*x) / exp(x);
    fprintf('r2 = %.10f\n\n', compute_fixed_point(g2, 1, f));
    
    fprintf('g3(x) = ((-7*x^2 + 3*x - e^x + 2) / -3.0)^(1/3) = x\n');
    g3 = @(x) ((-7*x^2 + 3*x - exp(x) + 2) / -3.0)^(1/3);
    fprintf('r3 = %0.10f\n\n', compute_fixed_point(g3, 2, f));
    
    fprintf('g4(x) = ln(3*x^3 - 7*x^2 + 3*x + 2) = x\n');
    g4 = @(x) log(3*x^3 - 7*x^2 + 3*x + 2);
    fprintf('r4 = %.10f\n', compute_fixed_point(g4, 6, f));
end


% Compute the fixed point of g(x).
function xc = compute_fixed_point(g, x0, f, tol)
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
