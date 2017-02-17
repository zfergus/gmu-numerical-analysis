% MATH 446: Project 01
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 01\nWritten by Zachary Ferguson\n\n');

    fprintf('Bisection Method:\n\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q1a
    fprintf('Q1a:\n\tf(x) = x^3 - 9\n');
    f = @(x) x^3 - 9;
    r = bisection_method(f, 2, 3);
    fprintf('\tr = %.10f\n', r);
    y = 3^(2. / 3.);
    print_errors(f, r, y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q1b
    fprintf('Q1b:\n\tf(x) = 3x^3 + x^2 - x - 5\n');
    f = @(x) 3 * x^3 + x^2 - x - 5;
    r = bisection_method(f, 1, 2);
    fprintf('\tr = %.10f\n', r);
    y = (1. / 9. * (-1 + (593 - 27 * (481^0.5))^(1.0 / 3.0) + ...
        (593 + 27 * (481^0.5))^(1.0 / 3.0)));
    print_errors(f, r, y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q1c
    fprintf('Q1c:\n\tf(x) = cos^2(x) - x + 6\n');
    f = @(x) cos(x) * cos(x) - x + 6;
    r = bisection_method(f, 6, 7);
    fprintf('\tr = %.10f\n', r);
    y = 6.7760923163195023262;
    print_errors(f, r, y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q3a
    fprintf('Q3a:\n\tf(x) = 2x^3 - 6x - 1\n');
    f = @(x) 2.0 * x.^3 - 6 * x - 1;

    x = linspace(-2, 2);
    y = f(x);
    figure;
    plot(x, y, x, 0*y);
    title('Q3a: f(x) = 2x^3 - 6x - 1');

    fprintf('\n\ta, b = -2, -1\n');
    r = bisection_method(f, -2, -1);
    fprintf('\tr1 = %.10f\n', r);
    y = -1.64178352745293;
    print_errors(f, r, y);

    fprintf('\n\ta, b = -1, 0\n');
    r = bisection_method(f, -1, 0);
    fprintf('\tr2 = %.10f\n', r);
    y = -0.168254401781027;
    print_errors(f, r, y);

    fprintf('\n\ta, b = 1, 2\n');
    r = bisection_method(f, 1, 2);
    fprintf('\tr3 = %.10f\n', r);
    y = 1.81003792923395;
    print_errors(f, r, y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q3b
    fprintf('Q3b:\n\tf(x) = e^(x-2) + x^3 - x\n');
    f = @(x) exp(x - 2) + x.^3 - x;

    x = linspace(-2, 2);
    y = f(x);
    figure;
    plot(x, y, x, 0*y);
    title('Q3b: f(x) = e^{x-2} + x^3 - x');

    fprintf('\n\ta, b = -2, -1\n');
    r = bisection_method(f, -2, -1);
    fprintf('\tr1 = %.10f\n', r);
    y = -1.0234821948582364944;
    print_errors(f, r, y);

    fprintf('\n\ta, b = -0.5, -0.5\n');
    r = bisection_method(f, -0.5, 0.5);
    fprintf('\tr2 = %.10f\n', r);
    y = 0.16382224325010849634;
    print_errors(f, r, y);

    fprintf('\n\ta, b = 0.5, 1.5\n');
    r = bisection_method(f, 0.5, 1.5);
    fprintf('\tr3 = %.10f\n', r);
    y = 0.78894138905554556637;
    print_errors(f, r, y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q3c
    fprintf('Q3c:\n\tf(x) = 1 + 5x - 6x^3 - e^(2x)\n');
    f = @(x) 1 + 5 * x - 6 * x.^3 - exp(2 * x);

    x = linspace(-1.5, 1.5);
    y = f(x);
    figure;
    plot(x, y, x, 0*y);
    title('Q3c: f(x) = 1 + 5x - 6x^3 - e^{2x}');

    fprintf('\n\ta, b = -1.5, -0.5\n');
    r = bisection_method(f, -1.5, -0.5);
    fprintf('\tr1 = %.10f\n', r);
    y = -0.81809373448119542124;
    print_errors(f, r, y);

    fprintf('\n\ta, b = -0.6, 0.4\n');
    r = bisection_method(f, -0.6, 0.4);
    fprintf('\tr2 = %.10f\n', r);
    y = 0.0;
    print_errors(f, r, y);

    fprintf('\n\ta, b = 0.5, 1.5\n');
    r = bisection_method(f, 0.5, 1.5);
    fprintf('\tr3 = %.10f\n', r);
    y = 0.50630828634622119599;
    print_errors(f, r, y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q4a
    fprintf('Q4a:\n\tf(x) = x^2 - A\n');
    fprintf('\tA = 2, (a, b) = (1, 2)\n');
    A = 2;
    f = @(x) x^2 - A;
    r = bisection_method(f, 1, 2);
    fprintf('\tr = %.10f\n', r);
    y = 2^0.5;
    print_errors(f, r, y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q4c
    fprintf('Q4b:\n\tA = 3, (a, b) = (1, 2)\n');
    A = 3;
    f = @(x) x^2 - A;
    r = bisection_method(f, 1, 2);
    fprintf('\tr = %.10f\n', r);
    y = 3^0.5;
    print_errors(f, r, y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Q4c
    fprintf('Q4c:\n\tA = 5, (a, b) = (2, 3)\n');
    A = 5;
    f = @(x) x^2 - A;
    r = bisection_method(f, 2, 3);
    fprintf('\tr = %.10f\n', r);
    y = 5^0.5;
    print_errors(f, r, y);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fixed Point Iteration
    fprintf('\nFixed Point Iteration:\n\n');

    fprintf('Q1a:\n\tg(x) = (2x+2)^(1/3) = x\n');
    fprintf('\txc = %.10f\n', fixed_point_iteration(...
        @(x)(2 * x + 2)^(1 / 3.), 2));

    fprintf('Q1b:\n\tg(x) = ln(7-x) = x\n');
    fprintf('\txc = %.10f\n', fixed_point_iteration(@(x) log(7 - x), 2));

    fprintf('Q1c:\n\tg(x) = ln(4-sin(x)) = x\n');
    fprintf('\txc = %.10f\n', ...
        fixed_point_iteration(@(x) log(4 - sin(x)), 2));

    A = 3.;
    g = @(x) (x + A / x) / 2.;
    fprintf('Q3a:\n\tg(x) = (x + 3 / x) / 2\n');
    fprintf('\tx0 = 2\n');
    fprintf('\txc = %.10f\n', fixed_point_iteration(g, 2));

    A = 5.;
    g = @(x) (x + A / x) / 2.;
    fprintf('Q3b:\n\tg(x) = (x + 5 / x) / 2\n');
    fprintf('\tx0 = 2\n');
    fprintf('\txc = %.10f\n', fixed_point_iteration(g, 2));

    fprintf('Q5:\n\tg(x) = cos^2(x)\n');
    xc = fixed_point_iteration(@(x) (cos(x))^2, 1, 1e-6);
    fprintf('\txc = %.10f\n', xc);
    fprintf('\n\td/dx g(x) = -2*cos(x)*sin(x)\n');
    fprintf('\t|d/dx g(xc)| = %.10f\n', abs(-2 * cos(xc) * sin(xc)));
    fprintf('\tTherefore g(x) is locally convergent to xc.\n');
end

% Print the forward and backward error of r.
function print_errors(f, r, y)
    fprintf('\tForward error: %.10f\n', abs(y - r));
    fprintf('\tBackward error: %.10f\n', abs(f(r)));
end
