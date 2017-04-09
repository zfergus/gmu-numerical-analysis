% MATH 446: Section 3.1 Question 1
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Section 3.1 Question 1\nWritten by Zachary Ferguson\n\n');

    % Data points for Section 3.1 Q1
    data = [1960, 3039585530; ...
            1970, 3707475887; ...
            1980, 4452584592; ...
            1990, 5281653820; ...
            2000, 6079603571];

    fprintf('=== Section 3.1 (Pg. 151) ===\n\nf(1980) = %d\n', data(3, 2));

    % Q1a
    fprintf('Q1a:\n');
    p1 = @(x) eval_newtdd(data, ...
        newtons_divided_differences([data(2, :); data(4, :)]), x);
    print_Q1(data(3, 2), p1, 'P_1');

    % Q1b
    fprintf('Q1b:\n');
    p2 = @(x) eval_newtdd(data, ...
        newtons_divided_differences([data(1:2, :); data(4, :)]), x);
    print_Q1(data(3, 2), p2, 'P_2');

    % Q1c
    fprintf('Q1c:\n');
    p3 = @(x) eval_newtdd(data, ...
        newtons_divided_differences([data(1:2, :); data(4:5, :);]), x);
    print_Q1(data(3, 2), p3, 'P_3');

    figure;
    plot(data(:, 1), data(:, 2), 'o');
    hold on;
    x = linspace(1950, 2010);
    plot(x, p1(x), '-r');
    plot(x, p2(x), '-g');
    plot(x, p3(x), '-b');
    hold off;

    fprintf('=== Section 3.2 (Pg. 157) ===\n\n');
end

function print_Q1(real_val, p, p_name)
    error = abs(real_val - p(1980));
    rel_err = error / real_val;
    fprintf('\t%s(1980) = %d\n\terror = %d\n\trelative error = %g\n', ...,
        p_name, int64(p(1980)), int64(error), rel_err);
end
