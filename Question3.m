
%constants
K = 1.3806e-23;
m = 0.26*9.1093e-31;

%Place electrons in Boundary
x = 200e-9*rand(1000,1);
y = 100e-9*rand(1000,1);


% Assign electron velocity based on Maxwell Boltzman Distribution
T = 300;
vth = sqrt(2*K*T/m);
std = sqrt(K*T/m);
Vx = normrnd(0,std,[1000,1]);
Vy = normrnd(0,std,[1000,1]);
V = sqrt(Vx.^2 + Vy.^2);
figure(1)
histogram(V);
xlabel('Velocity')
ylabel('Number of Electrons')
Tplot = zeros(1000,1); %calculated temperature for plotting

dt = 100e-9/vth/100;

%Initilizing varibales for calculating mean free path and mean time between
%collisions
numScat = zeros(1000,1); 
xScat = zeros(1000,1); % x distance between scatters
yScat = zeros(1000,1); % y distance between scatters
tScat = zeros(1000,1); % time between scatters
distScat = zeros(1000,1); % Total free distance travelled by each electron
totalfreetime = zeros(1000,1); % total time travlled between all scatters


for i =1:500
    xold = x;
    yold = y;
    
    % Define region boundaries and rules for interacting with boundaries
    xboundRight = x > 200e-9;
    xboundLeft = x < 0;
    ybound = (y > 100e-9) | (y <0);
    x(xboundRight) = x(xboundRight) - 200e-9;
    x(xboundLeft) = x(xboundLeft) + 200e-9;
    xold(xboundRight | xboundLeft) = x(xboundRight | xboundLeft);
    Vy(ybound) = -Vy(ybound);
    
    %Update Position
    x = x + Vx*dt;
    y = y + Vy*dt;
    
    % Determine Witchelectrons scatter and update velocity
    scatter = rand(1000,1) < (1 - exp(-dt/0.2e-12));
    Vx(scatter) = normrnd(0,std,size(Vx(scatter)));
    Vy(scatter) = normrnd(0,std,size(Vy(scatter)));
    
    %Determine distance and time between scatters for all electrons
    xScat(scatter) = x(scatter) - xScat(scatter);
    yScat(scatter) = y(scatter) - yScat(scatter);
    tScat(scatter) = dt*i - tScat(scatter);
    distScat(scatter) = distScat(scatter) + sqrt(xScat(scatter).^2 + yScat(scatter).^2);
    totalfreetime(scatter) = totalfreetime(scatter) + tScat(scatter);
    numScat(scatter) = numScat(scatter) + 1;

    xplot = transpose([xold(1:10) x(1:10)]);
    yplot = transpose([yold(1:10) y(1:10)]); 
    Tplot(i) = (1/(2*K))*mean(Vx.^2 + Vy.^2)*m;
    figure(2)
    subplot(2,1,1);
    plot(xplot,yplot)
    xlim([0 200e-9])
    ylim([0 100e-9])
    title('Electron Trajectory')
    xlabel('x')
    ylabel('y')
    hold on
    subplot(2,1,2)
    plot(Tplot(1:i))
    title('Temperature vs Time Step')
    xlabel('Number of Time Steps')
    ylabel('Temperature (K)')
    hold on
    drawnow
   
     
end

MFPelectron = distScat./numScat; %Mean Free Path of each electron
MFP = nanmean(MFPelectron) %overall Mean Free Path

tauElectron = totalfreetime./numScat; %Mean time between scatters for each electron
tau = nanmean(tauElectron) % Overall mean time between scatters