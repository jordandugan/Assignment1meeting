
K = 1.3806e-23;
m = 0.26*9.1093e-31;
T = 300;
vth = sqrt(2*K*T/m);
L = vth*2e-12;

x = 200e-9*rand(1000,1);
y = 100e-9*rand(1000,1);

theta = 2*pi*rand(1000,1);
Vx = vth.*cos(theta);
Vy = vth.*sin(theta);


dt = 100e-9/vth/100;


for i =1:500
    xold = x;
    yold = y;
    xboundRight = x > 200e-9;
    xboundLeft = x < 0;
    ybound = (y > 100e-9) | (y <0);
    x(xboundRight) = 0;
    x(xboundLeft) = 200e-9;
    xold(xboundRight|xboundLeft) = x(xboundRight|xboundLeft);
    Vy(ybound) = -Vy(ybound);
    x = x + Vx*dt;
    y = y + Vy*dt;
    
    xplot = transpose([xold(1:10) x(1:10)]);
    yplot = transpose([yold(1:10) y(1:10)]); 
    T(i) = (1/(2*K))*mean(Vx.^2 + Vy.^2)*m;
    subplot(2,1,1);
    plot(xplot,yplot)
    xlim([0 200e-9])
    ylim([0 100e-9])
    title('Electron Trajectory')
    xlabel('x')
    ylabel('y')
    hold on
    subplot(2,1,2)
    plot(T)
    title('Temperature vs Time Step')
    xlabel('Number of Time Steps')
    ylabel('Temperature (K)')
    hold on
    drawnow
    
end
