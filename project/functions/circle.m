%% Function which creates the round (circular) areas of the BS

function h = circle(x,y,r)

    hold on
    % theta angle (th)
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit);
    
end
