% Create a random set of coordinates in a circle.
% First define parameters that define the number of points and the circle.

function xy = punto_nel_cerchio(x0,y0,R,n)

    %R = 20;
    %x0 = 0; % Center of the circle in the x direction.
    %y0 = 0; % Center of the circle in the y direction.
    % Now create the set of points.
    t = 2*pi*rand(1,1);
    r = R*sqrt(rand(n,1));
    x = x0 + r.*cos(t);
    y = y0 + r.*sin(t);
    circle(x0,y0,R);
    xy=[x,y];
    ax=gca;
    figure(1)
    hold (ax,'on')
    % Now display our random set of points in a figure.
    plot(x,y, 'o', 'MarkerSize', 5)
    axis square;
    grid on;
    % Enlarge figure to full screen.
    %set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    %fontSize = 30;
    %xlabel('X', 'FontSize', fontSize);
    %ylabel('Y', 'FontSize', fontSize);
    %title('Random Locations Within a Circle', 'FontSize', fontSize);

end