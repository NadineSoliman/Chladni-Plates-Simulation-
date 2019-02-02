%Chladni Pattern Simulation For a Square Plate-- Project 3
%Introduction to Computer Smulations
%Nadine Soliman

%%% variables for equations %%%

% constants in equation
m=4;
n=3;
A =1;
B = 1;

%length of the square plate
l=1;
%velocity of propagation of the wave
v=2;

%spring constants
kx = m*pi/l;
ky = n*pi/l;

%angular velocity
w = v*sqrt(kx^2+ky^2);

% initialize x,y,z grid
x = -l:0.01:l;
y = -l:0.01:l;
z = zeros(size(x,2), size(y, 2));

% number of particles on the plate
N = 2500;

%%variables for plate with sand particles
%length of plate with sound x100 scale
L = 100;
kxL = m*pi/L;
kyL = n*pi/L;

%diameter of particles
diam = 0.5;

%drawing the gradient of the function
[uu,vv] = meshgrid(0:2:L,0:2:L);
Z = (A*sin(uu.*kxL).*sin(vv.*kyL)+B*sin(vv.*kxL).*sin(uu.*kyL));
[DX,DY] = gradient(Z,1,1);

% x,y coordinates of the sand particles
xx=zeros(1,N);
yy=zeros(1,N);

%space that could be occupied
space =L;

%initialize the positions of particles
for i=1:N
    xx(i) = rand*space;
    yy(i) = rand*space;
end

%axis limits setup
axisLimits2D = [-l l -l l];
axisLimits3D = [0 l 0 l -abs((A+abs(B))) abs((abs(A)+abs(B)))];

for t = 0:0.01:10000
    %make the plate vibrate
    for i =1:size(x,2)
        for j=1:size(y,2)
            z(i,j) = z(i,j)+(A*sin(kx*x(i))*sin(ky*y(j)) + B*sin(kx*y(j))*sin(ky*x(i)))*sin(w*t);
            z(i,j) = z(i,j)/(abs(A)+abs(B));
        end
    end
    
    %move the sand particles according to the gradient
    for i=1:N
        dx = kxL*A*cos(xx(i)*kxL)*sin(yy(i)*kyL)+B*kyL*cos(xx(i)*kyL)*sin(yy(i)*kyL);
        dy = kyL*A*sin(xx(i)*kxL)*cos(yy(i)*kyL)+B*kxL*sin(xx(i)*kyL)*cos(yy(i)*kxL);
        
        %move according to the negative of the gradient
        x_to_go = xx(i)- dx;
        y_to_go = yy(i)- dy;
        
%         %check for collisions
        for j=1:N
            if i~=j
                %if they hit make them bounce off eachother 
                if ((x_to_go -xx(j))^2 + (y_to_go - yy(j))^2)<(diam^2)
                    x_to_go =xx(i)+ (x_to_go -xx(j))*diam*0.25;
                    y_to_go =yy(i)+(y_to_go - yy(j))*diam*0.25;
                    break
                end
            end
        end
        %move the particles
        xx(i) = x_to_go;
        yy(i) = y_to_go;
    end
    
    %plotting
    %vibrating square palte
    subplot(2,2,1);
    s = surf(x,y,z);
    axis(axisLimits3D);
    
    %contour map
    subplot(2,2,2);
    contourf(x,y,z);
    axis(axisLimits2D);
    
    %predicted solution
    toPlot =  @(x,y) (A*sin(x.*kx).*sin(y.*ky) +B*sin(y.*kx).*sin(x.*ky));
    subplot(2,2,3);fimplicit(toPlot, [0 1 0 1]);
    
    %sand particles moving according to the vibrations
    subplot(2,2,4);
    %the flow field part
    quiver(uu,vv,DX.*sin(w*t),DY.*sin(w*t));
    hold on
    %the sand particles part
    scatter(xx, yy, 10);
    axis([0 L 0 L]);
    
    drawnow
    hold off
    drawnow
    axis manual
    axis equal
end
