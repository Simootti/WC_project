%% Cells initialization

% The 3 considered cells are circular and with radius R

% Setting the intersite distance, i.e. the distance between cells
% It's required to be: d=R*sqrt(3)
% So we first set the radius "R", for example R=2e3 and then "d"
d=R*sqrt(3);

% In order to limit the area in which the "ue" (users) could be positioned,
% we draw a square with border: ds>=d+2*R (or ds=2d)
ds=d+2*R;

% Axes origin is in the middle of the square
O=[ds,ds-R/5];

% Square with border: ds
x = [ds/2 ds/2 1.5*ds 1.5*ds ds/2]; 
y = [ds/2 1.5*ds 1.5*ds ds/2 ds/2];

shadow={};
% Implementation of 1:4 squares inside the main square in which
% users (devices) could be positionated
for i=1:4
    switch i
        case 1
                % Squares related to Shadowing
                shadow(i).x = [ds/2 ds/2 ds ds ds/2];
                shadow(i).y = [ds 1.5*ds 1.5*ds ds ds];
        case 2
                shadow(i).x = [ds ds 1.5*ds 1.5*ds ds];
                shadow(i).y = [ds 1.5*ds 1.5*ds ds ds];
        case 3
                shadow(i).x = [ds ds 1.5*ds 1.5*ds ds];
                shadow(i).y = [ds/2 ds ds ds/2 ds/2];
        case 4
                shadow(i).x = [ds/2 ds/2 ds ds ds/2];
                shadow(i).y = [ds/2 ds ds ds/2 ds/2];
    end
            
end

% We calculate the centers of the 3 cells thanks to the apotema's formula
% of the trianle exploiting the dimension of the side L=d : a = L/(2*sqrt(3))
L=d;
a = L/(2*sqrt(3));
% Y-axis of cells related to Bs 2 and 3 ( on the left on the right)
Y2=O(2)-a;
Y3=Y2;
% X-axis of cell related to Bs 2 
X2=O(1)-(L/2);
% X-axis of cell related to Bs 3
X3=O(1)+(L/2);
% X-axis of cell related to Bs 1
X1=O(1);
Y1=O(2)+(sqrt(L^2-(L/2)^2)-a);

% We draw the tringle which connects the Base Stations
tri = [X2 X1 X3 X2; Y2 Y1 Y3 Y2];

% We draw the circles corresponding to the 3 cells and the position
% of the 3 Base Stations
Bs={};

% Cella 1
Bs(1).pos = [X1 Y1]; 

% Left Cell, cell 2
Bs(2).pos = [X2 Y2]; 

% Rigth Cella, cell 3
Bs(3).pos = [X3,Y3];

figure
ax=gca;
hold (ax,'on')
%axis square
cella_superiore = circle(X1,Y1,R);
hold on
cella_sinistra = circle(X2,Y2,R);
cella_destra = circle(X3,Y3,R); 
plot(x,y,'r', 'LineWidth',3);
for i=1:4
    plot(shadow(i).x,shadow(i).y);
end
plot(O(1),O(2));
plot(tri(1,:), tri(2,:));
txt=["Bs1","Bs2","Bs3"];
for i=1:3
    plot(Bs(i).pos(1),Bs(i).pos(2),'*','DisplayName',txt(1,i))
end

hold off