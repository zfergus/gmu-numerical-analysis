% MATH 446: Project 02
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 02\nWritten by Zachary Ferguson\n\n');

    fprintf('f(x) = 3*x^3 - 7*x^2 + 3*x - e^x + 2 = 0\n\n');
    f = @(x) 3*x^3 - 7*x^2 + 3*x - exp(x) + 2;

    fprintf('g1(x) = (e^x - 2) / (3x^2 - 7x + 3) = x\n')
    g1 = @(x) (exp(x) - 2) / (3*x^2 - 7*x + 3);
    fprintf('r1 = %.10f\n\n', fixed_point_iteration(g1, -1, f));

    fprintf('g2(x) = (3*x^4 - 7*x^3 + 3*x^2 + 2*x) / e^x = x\n');
    g2 = @(x) (3*x^4 - 7*x^3 + 3*x^2 + 2*x) / exp(x);
    fprintf('r2 = %.10f\n\n', fixed_point_iteration(g2, 1, f));

    fprintf('g3(x) = ((-7*x^2 + 3*x - e^x + 2) / -3.0)^(1/3) = x\n');
    g3 = @(x) ((-7*x^2 + 3*x - exp(x) + 2) / -3.0)^(1/3);
    fprintf('r3 = %0.10f\n\n', fixed_point_iteration(g3, 2, f));

    fprintf('g4(x) = ln(3*x^3 - 7*x^2 + 3*x + 2) = x\n');
    g4 = @(x) log(3*x^3 - 7*x^2 + 3*x + 2);
    fprintf('r4 = %.10f\n', fixed_point_iteration(g4, 6, f));
end
