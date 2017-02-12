% Computes the roots of a function using the Newton's Method.
% Written by Zachary Ferguson


% Solve the project questions
function newtons_method
    % Q1a
    fprintf('Q1a:\n\tf(x) = x^3 - 2x - 2 = 0\n');
    f  = @(x) x^3 - 2*x - 2;
    fp = @(x) 3*x^2 - 2;
    x0 = 2;
    fprintf('\tx0 = %g\n', x0);
    fprintf('\tr = %.10f\n', compute_root(f, fp, x0));

    % Q1b
    fprintf('Q1b:\n\tf(x) = e^x + x - 7 = 0\n');
    f  = @(x) exp(x) + x - 7;
    fp = @(x) exp(x) + 1;
    x0 = 0;
    fprintf('\tx0 = %g\n', x0);
    fprintf('\tr = %.10f\n', compute_root(f, fp, x0));

    % Q1c
    fprintf('Q1c:\n\tf(x) = e^x + sin(x) - 4 = 0\n');
    f  = @(x) exp(x) + sin(x) - 4;
    fp = @(x) exp(x) + cos(x);
    x0 = 2;
    fprintf('\tx0 = %g\n', x0);
    fprintf('\tr = %.10f\n', compute_root(f, fp, x0));

    % Q3a
    fprintf('Q3a:\n\tf(x) = 27x^3 + 54x^2 + 36x + 8 = 0\n');
    f  = @(x) 27*x^3 + 54*x^2 + 36*x + 8;
    fp = @(x) 81*x^2 + 108*x + 36;
    x0 = 0.0;
    r = -2/3;
    fprintf('\tx0 = %g\n', x0);
    xc = compute_root(f, fp, x0, 1e-16);
    fprintf('\txc = %.16f\n', xc);
    fprintf('\tForward Error = |r - xc| = %.16f\n', abs(r-xc));
    fprintf('\tBackward Error = f(xc) = %.16f\n', f(xc));
    fprintf('\tmultiplicity of r is 3\n');
    xc = compute_root(f, fp, x0, 1e-16, 3);
    fprintf('\txc = %.16f\n', xc);
    fprintf('\tForward Error = |r - xc| = %.16f\n', abs(r-xc));
    fprintf('\tBackward Error = f(xc) = %.16f\n', f(xc));

    % Q9
    fprintf('Q9:\n\tf(x) = 14xe^(x-2) - 12e^(x-2) - 7x^3 + 20x^2 - 26x + 12\n');
    f  = @(x) 14*x*exp(x-2) - 12*exp(x-2) - 7*x^3 + 20*x^2 - 26*x + 12;
    fp = @(x) 14*x*exp(x-2) + 2*exp(x-2) - 21*x^2 + 40*x - 26;
    fpp = @(x) 14*x*exp(x-2) + 4*exp(x-2) - 42*x + 40
    x0 = 0;
    r = compute_root(f, fp, x0, 1e-9, 1, 2);
    fprintf('\tr1 = %.10f\n', r);
    fprintf('\tM = lim i->inf (e_(i+1)/(e_i)^2) = %.10f\n\n', ...
        abs(fpp(r)/(2*fp(r))))

    x0 = 3.0;
    r = compute_root(f, fp, x0, 1e-9, 1, 1);
    fprintf('\tr2 = %.10f\n', r);
    fprintf('\tmultiplicity of r2 is 3 -> ')
    fprintf('S = lim i->inf (e_(i+1)/e_i) = %.10f\n', 2/3)
end


% Compute the foot to f(x)
function r = compute_root(f, fp, x0, tol, m, print_ei)
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
    r = x;
end
