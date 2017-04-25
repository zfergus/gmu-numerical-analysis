% MATH 446: Project 11
% Written by Zachary Ferguson

function main()
    fprintf('MATH 446: Project 11\nWritten by Zachary Ferguson\n\n');

    fprintf('=== Section 3.2 (Pg. 157) ===\n\n');

    points = [1994, 67.052;
              1995, 68.008;
              1996, 69.803;
              1997, 72.024;
              1998, 73.400;
              1999, 72.063;
              2000, 74.669;
              2001, 74.487;
              2002, 74.065;
              2003, 76.777];

    coeffs = newtons_divided_differences(points);
    fprintf('--- Q3 ---\n\nPoints:\n');
    disp(points);
    p_str = newtdd_str(points, coeffs);
    fprintf('\nP(x) = %s\n', p_str);

    figure;
    x = linspace(1993, 2004, 200);
    y = eval_newtdd(points, coeffs, x);
    plot(x, y, '-m');
    hold on;

    fprintf('\nEstimate of oil production per day in 2010: %g\n', ...
        eval_newtdd(points, coeffs, 2010));

    fprintf('\nThis interpolation exhibits the Runge phenomenon.\n');
    fprintf('This interpolating polynomial is a bad model of the data\n');
    fprintf('because it does not model the data after the given points. \n');
    fprintf('This model does not transition smoothly from point to point.\n\n');

    fprintf('\n=== Section 3.4 (Pg. 178) ===\n\n');
    
    fprintf('--- Q13 ---\n\n');
    
    natural_coeffs = cubic_spline(points, 'natural');
    fprintf('Natural Spline:\n\n      a    |    b    |    c\n');
    disp(natural_coeffs);
    not_a_knot_coeffs = cubic_spline(points, 'not-a-knot');
    fprintf('\nNot-a-Knot Spline:\n\n      a     |     b    |    c\n');
    disp(not_a_knot_coeffs);
    parabola_coeffs = cubic_spline(points, 'parabola');
    fprintf('\nParaboloically Terminated Spline:\n\n      a    |    b    |    c\n');
    disp(parabola_coeffs);
    
    y = eval_cubic_spline(points, natural_coeffs, x);
    plot(x, y, '-r');
    y = eval_cubic_spline(points, not_a_knot_coeffs, x);
    plot(x, y, '-g');
    y = eval_cubic_spline(points, parabola_coeffs, x);
    plot(x, y, '-b');
    plot(points(:, 1), points(:, 2), 'ok')
    hold off;
    axis([1993 2004 60 95]);
    legend('Q_9(x)', 'Natural Spline', 'Not-a-Knot', 'Parabolically Terminated Spline');
    title('Total World Oil Production');
    xlabel('year');
    ylabel('bbl\day (x10^6)');
end
