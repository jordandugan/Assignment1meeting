K = 1.3806e-23;
m = 0.26*9.1093e-31;

x = 200e-9*rand(1000,1);
y = 100e-9*rand(1000,1);

yboundSpecular = true; %Specular reflection when treu, diffuse otherwize;
xboundSpecular = true;
boxSpecular = true;

inbox1 = x > 80e-9 & x < 120e-9 & y > 60e-9;
inbox2 = x > 80e-9 & x < 120e-9 & y < 40e-9;
x(inbox1) = x(inbox1) + ((rand() > 0.5)*2 - 1)*40e-9*rand(size(x(inbox1)));
x(inbox2) = x(inbox2) + ((rand() > 0.5)*2 - 1)*40e-9*rand(size(x(inbox2)));
y(inbox1) = y(inbox1) - 0.2*rand(size(y(inbox1)));
y(inbox2) = y(inbox2) + 0.2*rand(size(y(inbox2)));

T = 300;
vth = sqrt(2*K*T/m);
std = sqrt(K*T/m);
Vx = normrnd(0,std,[1000,1]);
Vy = normrnd(0,std,[1000,1]);
V = sqrt(Vx.^2 + Vy.^2);
Tplot = zeros(1000,1);
figure(1)
histogram(V);
dt = 0.5e-14;

xold = x;
yold = y;
for i =1:400
    %Defines the boundaries of the simulation as well as the boxes
    yboundTop = y > 100e-9;
    yboundBottom = y < 0;
    inbox1 = x >= 80e-9 & x <= 120e-9 & y >= 60e-9;
    inbox2 = x >= 80e-9 & x <= 120e-9 & y <= 40e-9;
    xboundRight = x > 200e-9;
    xboundLeft = x < 0;
    
    %Reflection off of xboundary
    if xboundSpecular
        Vx(xboundRight | xboundLeft) = - Vx(xboundRight | xboundLeft);
    else
        theta = pi*rand();
        Vx(xboundRight | xboundLeft) = V(xboundRight | xboundLeft)*cos(theta);
        Vy(xboundRight | xboundLeft) = V(xboundRight | xboundLeft)*sin(theta);
    end
    
    %Reflection off of y boundary
    if yboundSpecular
        Vy(yboundTop | yboundBottom) = -Vy(yboundTop | yboundBottom);
    else
        theta = pi*rand();
        Vy(yboundTop | yboundBottom) = V(yboundTop | yboundBottom)*cos(theta);
        Vx(yboundTop | yboundBottom) = V(yboundTop | yboundBottom)*sin(theta);
    end
    
    %Reflection off of box
    if boxSpecular
        %Reflection off of verticle face
        Vx(inbox1 & yold >= 60e-9) = -Vx(inbox1 & yold >= 60e-9);
        Vx(inbox2 & yold <= 40e-9) = -Vx(inbox2 & yold <= 40e-9);
        
        %Reflection off of Horizontal face
        Vy(inbox1 & yold <= 60e-9) = -Vy(inbox1 & yold <= 60e-9);
        Vy(inbox2 & yold >= 40e-9) = -Vy(inbox2 & yold >= 40e-9);
    else
        theta = pi*rand();
        
        %Reflection off of verticle face
        Vx(inbox1 & yold >= 60e-9) = V(inbox1 & yold >= 60e-9)*cos(theta);
        Vx(inbox2 & yold <= 40e-9) = V(inbox2 & yold <= 40e-9)*cos(theta);
        Vy(inbox1 & yold >= 60e-9) = V(inbox1 & yold >= 60e-9)*sin(theta);
        Vy(inbox2 & yold <= 40e-9) = V(inbox2 & yold <= 40e-9)*sin(theta);
        
        %Reflection off of Horizontal Face
        Vy(inbox1 & yold <= 60e-9) = V(inbox1 & yold <= 60e-9)*cos(theta);
        Vy(inbox2 & yold >= 40e-9) = V(inbox2 & yold >= 40e-9)*cos(theta);
        Vx(inbox1 & yold <= 60e-9) = V(inbox1 & yold <= 60e-9)*sin(theta);
        Vx(inbox2 & yold >= 40e-9) = V(inbox2 & yold >= 40e-9)*sin(theta);
    end
    
    %Make sure all electrons are in the proper boundary
    y(yboundTop) = 100e-9;
    y(yboundBottom) = 0;
    x(xboundRight) = 200e-9;
    x(xboundLeft) = 0;
    x(inbox1 & yold >= 60e-9 & x <= 100e-9) = 80e-9;
    x(inbox1 & yold >= 60e-9 & x > 100e-9) = 120e-9;
    x(inbox2 & yold <= 40e-9 & x <= 100e-9) =80e-9;
    x(inbox2 & yold <= 40e-9 & x >= 100e-9) =120e-9;
    y(inbox1 & yold <= 60e-9) = 60e-9;
    y(inbox2 & yold >= 60e-9) = 40e-9;
   

    xold = x;
    yold = y;
    x = x + Vx*dt;
    y = y + Vy*dt;
    
    %Scatter electrons randomly and reasign velocity using Maxwell Boltzman
    %distribution
    scatter = rand(1000,1) < (1 - exp(-dt/0.2e-12));
    Vx(scatter) = normrnd(0,std,size(Vx(scatter)));
    Vy(scatter) = normrnd(0,std,size(Vy(scatter)));
    
    
    xplot = transpose([xold(1:20) x(1:20)]);
    yplot = transpose([yold(1:20) y(1:20)]);
    figure(2)
    subplot(2,1,1)
    Tplot(i) = (1/(2*K))*mean(Vx.^2 + Vy.^2)*m;
    plot(xplot,yplot)
    xlim([0 200e-9])
    ylim([0 100e-9])
    hold on
    plot(1e-9*[80 80 120 120],1e-9*[200 60 60 200])
    plot(1e-9*[80 80 120 120],1e-9*[0 40 40 0])
    
    subplot(2,1,2)
    plot(Tplot(1:i))
    drawnow
    
    
   
     
end


;


figure(3);
hist3([x y],'CdataMode','auto'); 
view(2);
title('Electron Density');
xlabel('x (m)');
ylabel('y (m)');
title('Electron Density Heat Map');

temp_sum_x = zeros(20,10);
temp_sum_y = zeros(20,10);
temp_num = zeros(20,10);

for i=1:1000
 % Find which "bin" it belongs in:
 x1 = floor(x(i)/1e-8);
 y1 = floor(y(i)/1e-8);
 if(x1<=0)
 x1 = 1;
 end
 if(y1<=0)
 y1= 1;
 end
 if(y1>100)
     y1 = 100;
 end
 if(x1>200)
     x1=200;
 end
 % Add its velocity components to the cumulative count:
 temp_sum_y(x1,y1) = temp_sum_y(x1,y1) + Vy(i).^2;
 temp_sum_x(x1,y1) = temp_sum_x(x1,y1) + Vx(i).^2;
 temp_num(x1,y1) = temp_num(x1,y1) + 1;

end

temp = (temp_sum_x + temp_sum_y).*m./K./2./temp_num;
temp(isnan(temp)) = 0;
temp = transpose(temp);


figure(4)
surf(temp)
title('Temperature Map');
xlabel('x (nm)');
ylabel('y (nm)');
        
       
