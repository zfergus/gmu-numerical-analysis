% Calculates coefficents for a cubic spline
% Written by Zachary Ferguson

function coeffs = cubic_spline(points)
    % Calculates coefficients of cubic spline
    % Input:
    %   x,y vectors of data points plus two optional extra data v1, vn
    % Output:
    %   matrix of coefficients b1,c1,d1;b2,c2,d2;...

    n = size(points, 1);
    v1 = 0;
    vn = 0;


    deltas = points(2:n, :) - points(1:n-1, :);
    dx = deltas(:, 1);
    dy = deltas(:, 2);

    A = diag([dx(1:n-2); 0], -1) + ...
      diag([0; 2*(dx(1:n-2) + dx(2:n-1)); 0]) + ...
      diag([0; dx(2:n-1)], 1);

    r = [0; 3 * (dy(2:n-1)./dx(2:n-1) - dy(1:n-2)./dx(1:n-2)); 0]; % right-hand side

    % Set endpoint conditions
    % Use only one of following 5 pairs:
    A(1,1) = 1; % natural spline conditions
    A(n,n) = 1;
    %A(1,1)=2;r(1)=v1; % curvature-adj conditions
    %A(n,n)=2;r(n)=vn;
    %A(1,1:2)=[2*dx(1) dx(1)];r(1)=3*(dy(1)/dx(1)-v1); %clamped
    %A(n,n-1:n)=[dx(n-1) 2*dx(n-1)];r(n)=3*(vn-dy(n-1)/dx(n-1));
    %A(1,1:2)=[1 -1]; % parabol-term conditions, for n>=3
    %A(n,n-1:n)=[1 -1];
    %A(1,1:3)=[dx(2) -(dx(1)+dx(2)) dx(1)]; % not-a-knot, for n>=4
    %A(n,n-2:n)=[dx(n-1) -(dx(n-2)+dx(n-1)) dx(n-2)];

    coeffs = zeros(n,3);
    coeffs(:,2) = A\r; % solve for c coefficients
    coeffs(1:n-1,3) = (coeffs(2:n, 2) - coeffs(1:n-1, 2))./(3*dx);
    coeffs(1:n-1,1) = dy./dx - dx.*(2 * coeffs(1:n-1, 2) - coeffs(2:n, 2)) / 3;
    coeffs = coeffs(1:n-1, :);
end
