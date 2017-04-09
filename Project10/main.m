% MATH 446: Project 10
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 10\nWritten by Zachary Ferguson\n\n');

    fprintf('=== Section 3.1 (Pg. 151) ===\n\n');
    [p, cos1_str] = build_cos1();

    figure;
    x = linspace(-4*pi, 4*pi);
    y = p(x);
    plot(x, cos(x), '-bo');
    hold on;
    plot(x, y, '-rx');
    hold off;
    axis([-4*pi 4*pi -1.5 1.5]);
    legend('cos(x)', 'cos1(x)');
    title('Approximation of cos(x)');

    fprintf('Q4:\n')
    fprintf('\tcos1(x) = %s\n', cos1_str);
    fprintf('\tFundamental Domain of cos: [0, PI/2]\n');
    fprintf('\tSee Figure 1 for plot of cos1.\n');
    fprintf('\tforward error of cos1 = %g\n', norm(cos(x) - y, inf));
    fprintf('\tSee Figure 2 for plot of actual error of cos1(x).\n\n');

    figure;
    x = linspace(0, pi/2);
    y = abs(cos(x) - p(x));
    plot(x, y, '-r');
    title('Error of cos1(x)');


    fprintf('=== Section 3.2 (Pg. 157) ===\n\n');

    % Data points for Section 3.2 Q1
    data = [0.6, 1.433329; ...
            0.7, 1.632316; ...
            0.8, 1.896481; ...
            0.9, 2.247908; ...
            1.0, 2.718282];

    coeffs = newtons_divided_differences(data);

    p_str = newtdd_str(data, coeffs);
    fprintf('Q1a:\n\tP(x) = %s\n', p_str);

    p = @(x) eval_newtdd(data, coeffs, x);
    fprintf('Q1b:\n\tP(0.82) = %g\n\tP(0.98) = %g\n', p(0.82), p(0.98));

    f = @(x) exp(x.^2);
    upper_error1 = upper_limit_error(data, 312*exp(1), 0.82);
    upper_error2 = upper_limit_error(data, 312*exp(1), 0.98);
    fprintf('Q1c:\n');
    fprintf('\tupper limit of error @ x = 0.82: %g\n', upper_error1);
    fprintf('\tactual error @ x = 0.82: %g\n', abs(f(0.82) - p(0.82)));
    fprintf('\tupper limit of error @ x = 0.98: %g\n', upper_error2);
    fprintf('\tactual error @ x = 0.98: %g\n', abs(f(0.98) - p(0.98)));

    fprintf('Q1d:\n\tSee Figures 3 and 4 for plots of error.\n');
    figure;
    x1 = linspace(0.5, 1);
    y1 = abs(p(x1) - f(x1));
    plot(x1, y1, '-r');
    title('Actual Error for Range [0.5, 1]');

    figure;
    x2 = linspace(0, 2);
    y2 = abs(p(x2) - f(x2));
    plot(x2, y2, '-b');
    title('Actual Error for Range [0, 2]');
end

function ue = upper_limit_error(points, f_prime_c, x)
    n = size(points, 1);
    prod = 1;
    for i = 1:n
        prod = prod * (x - points(i, 1));
    end
    ue = abs(prod / factorial(n) * f_prime_c);
end
