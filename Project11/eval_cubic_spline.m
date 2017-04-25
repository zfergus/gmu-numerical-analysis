% Evaluates the given cubic spline at the x value.
% Written by Zachary Ferguson

function y = eval_cubic_spline(points, coeffs, x)
    % Evaluates the given cubic spline at the x value.
    % Input: x,y vectors of data points, coefficients of spline, x value(s)
    % Output: y values for given x values
    y = [];

    np = size(points, 1);
    nx = size(x, 2);
    for i = 1:nx
        if x(i) < points(1, 1)
            px = points(1, 1);
            j = 1;
        elseif x(i) >= points(end, 1)
            px = points(np-1, 1);
            j = np - 1;
        else
            for j = 2:np
                if x(i) < points(j, 1)
                    px = points(j-1, 1);
                    j = j - 1;
                    break;
                end
            end
        end
        dx = x(i) - px;
        yi = coeffs(j,3)*dx; % evaluate using nested multiplication
        yi = (yi+coeffs(j,2)).*dx;
        yi = (yi+coeffs(j,1)).*dx + points(j, 2);
        y  = [y; yi];
    end
end
