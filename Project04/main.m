% MATH 446: Project 04
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 04\nWritten by Zachary Ferguson\n\n');

    init_guess_str = 'x0 = 1, x1 = 2';
    x0 = 1;
    x1 = 2;

    f1_str = 'f(x) = x^3 - 2x - 2 = 0';
    f1 = @(x) x^3 - 2*x - 2;

    f2_str = 'f(x) = e^x + x - 7 = 0';
    f2 = @(x) exp(x) + x - 7;

    f3_str = 'f(x) = e^x + sin(x) - 4 = 0';
    f3 = @(x) exp(x) + sin(x) - 4;

    fprintf('Secant Method:\n\n');

    % Q1a
    fprintf('Q1a:\n\t%s\n', f1_str);
    fprintf('\t%s\n', init_guess_str);
    fprintf('\txc = %.10f\n', secant_method(f1, x0, x1));

    % Q1b
    fprintf('Q1b:\n\t%s\n', f2_str);
    fprintf('\t%s\n', init_guess_str);
    fprintf('\txc = %.10f\n', secant_method(f2, x0, x1));

    % Q1a
    fprintf('Q1c:\n\t%s\n', f3_str);
    fprintf('\t%s\n', init_guess_str);
    fprintf('\txc = %.10f\n', secant_method(f3, x0, x1));

    fprintf('\nMethod of False Position:\n\n');

    % Q2a
    fprintf('Q2a:\n\t%s\n', f1_str);
    fprintf('\t%s\n', init_guess_str);
    fprintf('\txc = %.10f\n', method_of_false_position(f1, x0, x1));

    % Q2b
    fprintf('Q2b:\n\t%s\n', f2_str);
    fprintf('\t%s\n', init_guess_str);
    fprintf('\txc = %.10f\n', method_of_false_position(f2, x0, x1));

    % Q2c
    fprintf('Q2c:\n\t%s\n', f3_str);
    fprintf('\t%s\n', init_guess_str);
    fprintf('\txc = %.10f\n', method_of_false_position(f3, x0, x1));

    fprintf('\nInverse Quadratic Interpolation:\n\n');
    
    x2 = 0;
    
    % Q2a
    fprintf('Q3a:\n\t%s\n', f1_str);
    fprintf('\t%s, x2 = %d\n', init_guess_str, x2);
    fprintf('\txc = %.10f\n', ...
        inverse_quadratic_interpolation(f1, x0, x1, x2));

    % Q2b
    fprintf('Q3b:\n\t%s\n', f2_str);
    fprintf('\t%s, x2 = %d\n', init_guess_str, x2);
    fprintf('\txc = %.10f\n', ...
        inverse_quadratic_interpolation(f2, x0, x1, x2));

    % Q2c
    fprintf('Q3c:\n\t%s\n', f3_str);
    fprintf('\t%s, x2 = %d\n', init_guess_str, x2);
    fprintf('\txc = %.10f\n', ...
        inverse_quadratic_interpolation(f3, x0, x1, x2));
end
