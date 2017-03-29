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

    titles = ['A', 'B'];
    figures = [];
    for i = 1:4
        figures = [figures figure];
        title(sprintf('Part %s (Solution %d): Comparison of Methods', ...
            titles(ceil(i / 2)), mod(i-1, 2) + 1));
        xlabel('Iteration Step');
        ylabel('Forward error of x_k');
        set(gca, 'YScale', 'log');
        hold on;
    end

    tol = 1e-10;

    % Q5a
    fprintf('Q5a:\n\tIntersection Point 1:\n');
    x0_a1 = zeros(3, 1);
    fprintf('\tx0 = \n');
    disp(x0_a1);
    xc = multivariate_newtons_method(f_a, df_a, x0_a1, tol, figures(1));
    fprintf('\txc = \n');
    disp(xc);

    fprintf('\tIntersection Point 2:\n');
    x0_a2 = 2*ones(3, 1);
    fprintf('\tx0 = \n');
    disp(x0_a2);
    xc = multivariate_newtons_method(f_a, df_a, x0_a2, tol, figures(2));
    fprintf('\txc = \n');
    disp(xc);

    % Q5b
    fprintf('Q5b:\n\tIntersection Point 1:\n');
    x0_b1 = -2*ones(3, 1);
    fprintf('\tx0 = \n');
    disp(x0_b1);
    xc = multivariate_newtons_method(f_b, df_b, x0_b1, tol, figures(3));
    fprintf('\txc = \n');
    disp(xc);

    fprintf('\tIntersection Point 2:\n');
    x0_b2 = 2*ones(3, 1);
    fprintf('\tx0 = \n');
    disp(x0_b2);
    xc = multivariate_newtons_method(f_b, df_b, x0_b2, tol, figures(4));
    fprintf('\txc = \n');
    disp(xc);

    % Q9a
    fprintf('Q9a:\n\tIntersection Point 1:\n');
    xc = broydens_method_1(f_a, eye(3), x0_a1, tol, figures(1));
    fprintf('\txc = \n');
    disp(xc);

    fprintf('\tIntersection Point 2:\n');
    xc = broydens_method_1(f_a, eye(3), x0_a2, tol, figures(2));
    fprintf('\txc = \n');
    disp(xc);

    % Q9b
    fprintf('Q9b:\n\tIntersection Point 1:\n');
    xc = broydens_method_1(f_b, eye(3), x0_b1, tol, figures(3));
    fprintf('\txc = \n');
    disp(xc);

    fprintf('\tIntersection Point 2:\n');
    xc = broydens_method_1(f_b, eye(3), x0_b2, tol, figures(4));
    fprintf('\txc = \n');
    disp(xc);

    % Q11a
    fprintf('Q11a:\n\tIntersection Point 1:\n');
    xc = broydens_method_2(f_a, eye(3), x0_a1, tol, figures(1));
    fprintf('\txc = \n');
    disp(xc);

    fprintf('\tIntersection Point 2:\n');
    xc = broydens_method_2(f_a, eye(3), x0_a2, tol, figures(2));
    fprintf('\txc = \n');
    disp(xc);

    % Q11b
    fprintf('Q11b:\n\tIntersection Point 1:\n');
    xc = broydens_method_2(f_b, eye(3), x0_b1, tol, figures(3));
    fprintf('\txc = \n');
    disp(xc);

    fprintf('\tIntersection Point 2:\n');
    xc = broydens_method_2(f_b, eye(3), x0_b2, tol, figures(4));
    fprintf('\txc = \n');
    disp(xc);

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
