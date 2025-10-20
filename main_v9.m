% Computes spectral energy density versus wavenumber
% Published with "A procedure for producing energy spectra in riverine settings from gridded repeat vessel-based ADCP transects"
% under AGPL License
% 
% uVel, vVel, wVel = streamwise, transverse, vertical velocity matrices
% Ud, Vd, Wd = streamwise, transverse, vertical velocity vectors at nearest depth to surface

% Name of VMT composite transect file
filename='TRTS20230817.mat';

close all

load(filename)

hFig1 = figure;
uVel = V.uSmooth/100;
surface(V.mcsDist, V.mcsDepth, uVel,'AlignVertexCenters','on','FaceColor','interp','EdgeColor','none')
xlabel('Distance [m]')
ylabel('Depth [m]')
ax = gca;
ax.YDir = 'reverse'; 
pbaspect([1.5,1,1])
set(hFig1, 'Position', [0 0 1500 1100])
colormap(jet)
colorbar
title('Streamwise Velocity [m/s]')
set(gca,'FontSize',25);
exportgraphics(gcf,[strtok(filename,'.'),'_Streamwise Velocity.jpg']);

hFig2 = figure;
vVel = V.vSmooth/100;
surface(V.mcsDist, V.mcsDepth,vVel,'AlignVertexCenters','on','FaceColor','interp','EdgeColor','none')
xlabel('Distance [m]')
ylabel('Depth [m]')
ax = gca;
ax.YDir = 'reverse'; 
pbaspect([1.5,1,1])
set(hFig2, 'Position', [0 0 1500 1100])
colormap(jet)
colorbar
clim([-0.5 0.5])
title('Transverse Velocity [m/s]')
set(gca,'FontSize',25);
exportgraphics(gcf,[strtok(filename,'.'),'_Transverse Velocity.jpg']);

hFig3 = figure;
wVel = V.wSmooth/100;
surface(V.mcsDist, V.mcsDepth,wVel,'AlignVertexCenters','on','FaceColor','interp','EdgeColor','none')
xlabel('Distance [m]')
ylabel('Depth [m]')
ax = gca;
ax.YDir = 'reverse'; 
pbaspect([1.5,1,1])
set(hFig3, 'Position', [0 0 1500 1100])
colormap(jet)
colorbar
clim([-0.12 0.12])
title('Vertical Velocity [m/s]')
set(gca,'FontSize',25);
exportgraphics(gcf,[strtok(filename,'.'),'_Vertical Velocity.jpg']);

% bin selection for depth
wave_number_bin = 1;
wave_number_depth = V.mcsDepth(wave_number_bin,1);
fprintf('Depth for wavenumber computation %f\n',wave_number_depth)
Ud = uVel(wave_number_bin,:);
Vd = vVel(wave_number_bin,:);
Wd = wVel(wave_number_bin,:);

% Plot the velocity traces
dx = A.hgns;      % spatial discretization (typically 1.0 meter), this will determine the spatial frequency
dz = A.vgns;      % spatial discretization (typically 0.4 meter)
Fs = 1/dx;   % spatial sampling frequency
L=length(Ud); % Number of samples
t=(0:L-1)/(Fs); % Time vector
hFig1 = figure; 
set(hFig1, 'Position', [0 0 1500 1100])
plot(t,Ud,'LineWidth', 2); % Plot the velocity trace
hold on
plot(t,Vd,'LineWidth', 2);
plot(t,Wd,'LineWidth', 2);
title('Velocity Trace');
ylabel('Velocity [m/s]');
xlabel('time (s)');
legend('X','Y','Z')
set(gca,'FontSize',25);
pbaspect([1.5,1,1])
exportgraphics(gcf,[strtok(filename,'.'),'_Velocity Trace.jpg']);     

%% Plot the energy density spectrum (Frequency)
Ud_detrended = detrend(Ud,3);
n=L/2; % FFT will yield half the number of unique points
aFreq=Fs*(1:n)/n; % Frequency array (half the length of signal)
PowerU = abs(fft(Ud_detrended)).^2/L; 
% Plot of Power Density in [m^2/s^2] = [J/kg]
% Dividing Power density over Fs [1/m] we get Energy density 
% Plot of Energy Density in [m^3/s^2] = [J⋅m/kg]
EnergyU=PowerU/Fs;
hFig4 = figure; 
set(hFig4, 'Position', [0 0 1500 1100])
loglog(aFreq(2:end),EnergyU(2:L/2),'LineWidth', 2)
hold on
xline(0.1,':','Color','red','LineWidth', 2);
% plot for k^(-5/3)
range = [11e-2 1];
loglog(range, 0.007*range.^(-5/3), 'g', 'LineWidth', 4);
pbaspect([1.5,1,1])
xlim([1e-3 1])
ylim([1e-6 2])
title('Energy Density Spectrum')
xlabel ('Wavenumber [1/m]')
ylabel ('Energy Density [J⋅m/kg]') 
legend('Streamwise','Cutoff for E(k)','k^-^5^/^3')
set(gca,'FontSize',25);
exportgraphics(gcf,[strtok(filename,'.'),'_Energy Density Spectrum.jpg']);     

% Calculate the integral length scale [in seconds]
% [m^2/s^2].[m/s].[s^2/m^2]=[s]
uprimeU=std(Ud);
UbarU=mean(Ud);
%%
% the elements 2 to 10 are chosen as they asymptotes to a fixed number at the lowest wavenumbers
Ef0U=mean(EnergyU(2:8))    
IntegralLengthScaleU=(Ef0U*UbarU)/(4*uprimeU.^2)
% multiplying integral length scale [s] by mean streamwise velocity [m/s]
% gives us the eddy size in meters
disp(['Streamwise Integral Length Scale [m]: ',num2str(IntegralLengthScaleU*UbarU)])
