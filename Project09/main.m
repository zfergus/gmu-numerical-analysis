% MATH 446: Project 09
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 09\nWritten by Zachary Ferguson\n\n');

    % Function and Jacobian used in part a:
    [f_a, df_a] = build_funtion_and_jacobian([1,1,0], 1, [1,0,1], 1, ...
        [0,1,1], 1);
    % Function and Jacobian used in part b:
    [f_b, df_b] = build_funtion_and_jacobian([1,-2,0], 5, [-2,2,-1], 5, ...
        [4,-2,3], 5);

    x0s_a = [  zeros(3, 1), 2*ones(3, 1)];
    x0s_b = [-2*ones(3, 1), 2*ones(3, 1)];

    x_stars_a = [1/3 * ones(3, 1), ones(3, 1)];
    x_stars_b = [(1:3)', [17/9; 22/9; 19/9]];

    titles = ['A', 'B'];
    figures = [];
    for i = 1:4
        figures = [figures figure];
        title(sprintf('Part %s (Solution %d): Comparison of Methods', ...
            titles(ceil(i / 2)), mod(i-1, 2) + 1));
        xlabel('Iteration Step');
        ylabel('Backwards error of x_k');
        set(gca, 'YScale', 'log');
        hold on;
    end

    tol = 1e-12;

    % Q5a
    run_method('5a', f_a, df_a, x0s_a, x_stars_a, ...
        'multivariate_newtons_method', tol, figures(1:2));

    % Q5b
    run_method('5b', f_b, df_b, x0s_b, x_stars_b, ...
        'multivariate_newtons_method', tol, figures(3:4));

    % Q9a
    run_method('9a', f_a, eye(3), x0s_a, x_stars_a, 'broydens_method_1', ...
        tol, figures(1:2));

    % Q9b
    run_method('9b', f_b, eye(3), x0s_b, x_stars_b, 'broydens_method_1', ...
        tol, figures(3:4));

    % Q11a
    run_method('11a', f_a, eye(3), x0s_a, x_stars_a, 'broydens_method_2', ...
        tol, figures(1:2));

    % Q11b
    run_method('11b', f_b, eye(3), x0s_b, x_stars_b, 'broydens_method_2', ...
        tol, figures(3:4));

    for i = 1:4
        figure(figures(i));
        legend(['Multivariate Newton''' 's Method'], ...
            ['Broyden''' 's Method I'], ...
            ['Broyden''' 's Method II']);
        hold off;
    end
end

function [f, df] = build_funtion_and_jacobian(c1, r1, c2, r2, c3, r3)
    % Build a function f(x) for the intersection of three circles.
    % Helper function for building f and df in Q5.
    % Function
    f = @(x) [(x(1)-c1(1))^2 + (x(2)-c1(2))^2 + (x(3)-c1(3))^2 - r1^2; ...
              (x(1)-c2(1))^2 + (x(2)-c2(2))^2 + (x(3)-c2(3))^2 - r2^2; ...
              (x(1)-c3(1))^2 + (x(2)-c3(2))^2 + (x(3)-c3(3))^2 - r3^2];
    % Jacobian
    df = @(x) [2*(x(1)-c1(1)) 2*(x(2)-c1(2)) 2*(x(3)-c1(3)); ...
               2*(x(1)-c2(1)) 2*(x(2)-c2(2)) 2*(x(3)-c2(3)); ...
               2*(x(1)-c3(1)) 2*(x(2)-c3(2)) 2*(x(3)-c3(3))];
end

function run_method(q_num, f, df, x0s, x_stars, method, tol, figures)
    fprintf('=== Q%s: ===\n\n', q_num);
    for i = 1:size(x0s, 2)
        fprintf('--- Solution Point %d: ---\nx0 = \n', i);
        disp(x0s(:, i));
        xc = feval(method, f, df, x0s(:, i), tol, figures(i));
        fprintf('xc = \n');
        disp(xc);

        fprintf('Backwards Error = %g\n', norm(f(xc), inf));
        fprintf('Forwards Error = %g\n\n', norm(xc - x_stars(:, i), inf));
    end
end
