% Idea of this script is to fit some function to the ZLP of a spectrum,
% then remove a scaled fraction of this function from the ZLP.
%% Step 0: Load a Linescan into Matlab, and call it 'Linescan'.
clear EN1;
clear EN2;
clear fitwin1;
clear fitwin2;
% clear EN;
clear fitwin;
clear fitEN;
clear SD;
clear RChiSq;
clear Residual;
clear fitfunction;
clear fitfunctionwin;

%% Step 1: Define Energy Scale and Size of Spectrum. 
% Define variable name for input and string name for tiff file to save data
% to.  
Input = Linescan_N; LinescanRaw=Linescan_N;
Output = Input;
String = 'LinescanNicanite_PearsonVII_Backsub';

% % Define the energy scale for the linescan. Comment/Uncomment if you want
% % to include or skip this step. 
% startEn = -0.223; % Set the starting energy for the spectrum
% dispersion = 0.0005; % 
% sizeen=size(Input,2); % Calculate size of spectrum
% sizex=size(Input,1); % Calculate size of linescan
% EN=(startEn:dispersion:startEn+dispersion*(sizeen-1)); % Define the energy axis
% Channel0 = round(1+(-startEn/dispersion)); % IMPORTANT - This calculates the channel for Zero Energy Loss.

Scalefactor = 1; % Define scale factor for background subtraction. Normally set to 1, but can be set to less than 1 if you want to remove partial background before trying another function. 

% Define the sizes of the windows for background fitting
% Windows = [Channel0+25,Channel0+40,Channel0+230,Channel0+290];
% Windows = [Channel0+50,Channel0+85,Channel0+220,Channel0+280];
% Windows = [Channel0+13,Channel0+18,Channel0+90,Channel0+120];
Windows = [Channel0+14,Channel0+18,Channel0+90,Channel0+120];
EN1 = EN(Windows(1):Windows(2));
EN2 = EN(Windows(3):Windows(4));

% Choose spectra to view in plots
PlotBuffer = 20; % Choose how far (in pixels) beyond your windows you would like to plot data. 
Viewspec1 = 29;
Viewspec2 = 36;

% Specify fit parameters
E0 = 0.0; % Min position for centre of ZLP
E0Var =   0.05; % Max position for centre of ZLP
FWHM = 0.02; % FWHM of Pearson VII function Formula
FWHMVar = 0.02; % Allowed variation in FWHM.
Powmin = 1; % Minimum power of Pearson VII function
Powmax = 3; % Maximum power of Pearson VII function
Powguess = 2; % Initial guess for power of Pearson VII function

options = optimoptions('lsqcurvefit','Display','none', 'FiniteDifferenceType','central', 'MaxFunctionEvaluations',1800,'MaxIterations', 10000);


%% Step 2: Fit Function, Background Subtract, Define Fit Quality
RChiSq = zeros(1,sizex);
RSq = zeros(1,sizex);

for m = 1:sizex

% % Comment/Uncomment this to switch on/off a 1 pixel median filter    
% sizeen=int16(sizeen);
% Output(m,1)= median(Output(m,1:2));
% for n = 2:sizeen-1
%     Output(m,n)= median([Output(m,n),Output(m,n-1),Output(m,n+1)]);% Pick out the mth spectrum for the curve fitting. 
% end
% Output(m,sizeen)= median(Output(m,sizeen-1:sizeen));
%     
% Comment/Uncomment this to switch to no oversampling    
FitSpec = (Output(m,:));

% % Comment/Uncomment this to switch to oversample by 1 pixel    
% if m == 1    
%     FitSpec = (Linescan_TailSubtract(m,:)+Linescan_TailSubtract(m+1,:))/2; % Pick out the mth spectrum for the curve fitting. 
% elseif m == sizex
%     FitSpec = (Linescan_TailSubtract(m,:)+Linescan_TailSubtract(m-1,:))/2; % Pick out the mth spectrum for the curve fitting.
% else
%     FitSpec = (Linescan_TailSubtract(m,:)+Linescan_TailSubtract(m+1,:)+Linescan_TailSubtract(m-1,:))/3; % Pick out the mth spectrum for the curve fitting. 
% end

% Define two windows to fit the ZLP. This allows you to leave the phonons
% out of the fitting window in case they skew the fit. 

fitwin1=FitSpec(Windows(1):Windows(2));
fitwin2=FitSpec(Windows(3):Windows(4));
fitwin = cat(2,fitwin1,fitwin2); %Concatenates the values from the two windows. 

Rawwin1=LinescanRaw(m,Windows(1):Windows(2));
Rawwin2=LinescanRaw(m,Windows(3):Windows(4));
Rawwin = cat(2,Rawwin1,Rawwin2); %Concatenates the values from the two windows. 

fitEN = cat(2,EN1,EN2);

%Fit Pearson VII (basically Lorentzian to a power) to the peaks.
func= @(C,x)PearsonVII(x,C(1),C(2),C(3),C(4))+C(5);
C1 = [max(FitSpec),          E0,           FWHM/2, Powguess, median(fitwin2)]; % Initial guesses for the fit parameters
lb = [max(FitSpec)/2,  E0-E0Var, (FWHM-FWHMVar)/2,   Powmin, min(FitSpec)]; % Lower bounds for the fit parameters
ub = [max(FitSpec)*40, E0+E0Var, (FWHM+FWHMVar)/2,   Powmax, max(fitwin2)]; % Upper bounds for the fit parameters 
%FR = linspace(fitEN(1),fitEN(end));
%figure(1);
ft = lsqcurvefit(func,C1,fitEN,fitwin, lb, ub, options); % Fit Curve to data.
% plot(fitEN,fitwin,'ro',FR,func(ft,FR),'b-');

% Make a measurement of the quality of the fit function within the window
fitfunctionwin(1,:) = PearsonVII(fitEN(1,:),ft(1),ft(2),ft(3),ft(4))+ft(5);
% SD(1,:) = abs(sqrt(Rawwin(1,:)));
% Residual(1,:) = ((fitwin(1,:)-fitfunctionwin(1,:))./SD(1,:)).^2;
% Estimate standard deviation from noise in low energy channels
SD = std(Input(m,1:10));
Residual=zeros(1,size(fitwin,2));
Res=zeros(1,size(fitwin,2));
Residual(1,:) = fitwin(1,:)-fitfunctionwin(1,:);
Res(1,:) = (Residual(1,:)./SD).^2;
RChiSq(1,m) = sum(Res)/size(fitEN,2);
%RChiSq(1,m)

% Make a measurement of the quality of the fit function within the window
% using an R squared metric.
R2 = RSquare1D(fitwin(1,:), Residual);
RSq(1,m)=R2;

% Subtract scaled fitted function from the data.  
fitfunction(1,:) = PearsonVII(EN(1,:),ft(1),ft(2),ft(3),ft(4))+ft(5);
Output(m,:) = Input(m,:) - Scalefactor*fitfunction;
Output(m,1:Windows(1)-1)=0;
% LinescanZLPPosition = ft2(2); % Spits out position of peak.

end
%% Step 2: Plot subtracted data.
close all;

[MRChiSq] = median(RChiSq);
[MaxRChiSq, MaxInd] = max(RChiSq);
[MinRChiSq, MinInd] = min(RChiSq);

figure('Name','Spectra_and_Backgrounds','units','normalized','outerposition',[0 0 1 1]);
set(gcf, 'color', [1 1 1]);
subplot(2,2,1)
FitSpec = (Input(Viewspec1,:));
fitwin1=FitSpec(Windows(1):Windows(2));
fitwin2=FitSpec(Windows(3):Windows(4));
fitwin = cat(2,fitwin1,fitwin2); %Concatenates the values from the two windows. 
fitEN = cat(2,EN1,EN2);
%Fit Gaussian plus Lorentzian to the peaks.
%Fit Gaussian plus Lorentzian to the peaks.
func= @(C,x)PearsonVII(x,C(1),C(2),C(3),C(4))+C(5);
C1 = [max(FitSpec),          E0,           FWHM/2, Powguess, median(fitwin2)]; % Initial guesses for the fit parameters
lb = [max(FitSpec)/2,  E0-E0Var, (FWHM-FWHMVar)/2,   Powmin, min(FitSpec)]; % Lower bounds for the fit parameters
ub = [max(FitSpec)*40, E0+E0Var, (FWHM+FWHMVar)/2,   Powmax, max(fitwin2)]; % Upper bounds for the fit parameters 
FR = linspace(fitEN(1),fitEN(end));
ft = lsqcurvefit(func,C1,fitEN,fitwin, lb, ub, options);
hold on;
plot(EN,Input(Viewspec1,:),'ro', FR,func(ft,FR),'b-');
fitfunc1P7=func(ft,FR); % Outputs the fit as an array
fitfunc1P7=fitfunc1P7';
title(['Spectrum ' num2str(Viewspec1) ': ChiSq ' num2str(RChiSq(Viewspec1)) ', Median ChiSq ' num2str(MRChiSq)]);
xlim([EN(Windows(1)-PlotBuffer) EN(Windows(4)+PlotBuffer)]);
ylim([0 fitwin(1)]);
ax = gca; ax.YAxis.Exponent = 0;
x1 = EN(Windows(1));
x2 = EN(Windows(2));
y2 = [min(Input(Viewspec1,:)) max(Input(Viewspec1,:))];
plot([x1 x1],y2, 'b-')
plot([x2 x2],y2, 'b-')    
x3 = EN(Windows(3));
x4 = EN(Windows(4));
y4 = [min(Input(Viewspec1,:)) max(Input(Viewspec1,:))];
plot([x3 x3],y4, 'g-')
plot([x4 x4],y4, 'g-') 
box on;
hold off;

subplot(2,2,2)
FitSpec = (Input(Viewspec2,:));
fitwin1=FitSpec(Windows(1):Windows(2));
fitwin2=FitSpec(Windows(3):Windows(4));
fitwin = cat(2,fitwin1,fitwin2); %Concatenates the values from the two windows. 
fitEN = cat(2,EN1,EN2);
%Fit Gaussian plus Lorentzian to the peaks.
%Fit Gaussian plus Lorentzian to the peaks.
func= @(C,x)PearsonVII(x,C(1),C(2),C(3),C(4))+C(5);
C1 = [max(FitSpec),          E0,           FWHM/2, Powguess, median(fitwin2)]; % Initial guesses for the fit parameters
lb = [max(FitSpec)/2,  E0-E0Var, (FWHM-FWHMVar)/2,   Powmin, min(FitSpec)]; % Lower bounds for the fit parameters
ub = [max(FitSpec)*40, E0+E0Var, (FWHM+FWHMVar)/2,   Powmax, max(fitwin2)]; % Upper bounds for the fit parameters 
FR = linspace(fitEN(1),fitEN(end));
ft = lsqcurvefit(func,C1,fitEN,fitwin, lb, ub, options); % Fit Curve to data.
hold on;
plot(EN,Input(Viewspec2,:),'ro', FR,func(ft,FR),'b-');
fitfunc2P7=func(ft,FR); % Outputs the fit as an array
fitfunc2P7=fitfunc2P7';
title(['Spectrum ' num2str(Viewspec2) ': ChiSq ' num2str(RChiSq(Viewspec2)) ', Median ChiSq ' num2str(MRChiSq)]);
xlim([EN(Windows(1)-PlotBuffer) EN(Windows(4)+PlotBuffer)]);
ylim([min(fitwin) max(fitwin)]);
ax = gca; ax.YAxis.Exponent = 0;
x1 = EN(Windows(1));
x2 = EN(Windows(2));
y2 = [min(Input(Viewspec2,:)) max(Input(Viewspec2,:))];
plot([x1 x1],y2, 'b-')
plot([x2 x2],y2, 'b-')    
x3 = EN(Windows(3));
x4 = EN(Windows(4));
y4 = [min(Input(Viewspec2,:)) max(Input(Viewspec2,:))];
plot([x3 x3],y4, 'g-')
plot([x4 x4],y4, 'g-') 
box on;
hold off;

subplot(2,2,3)
plot(EN,Output(Viewspec1,:),'r-');
title(['Spectrum ' num2str(Viewspec1) ':, ChiSq ' num2str(RChiSq(Viewspec1)) ', Median ChiSq ' num2str(MRChiSq)]);
xlim([EN(Windows(1)-PlotBuffer) EN(Windows(4)+PlotBuffer)]);
ylim([1.2*min(Output(Viewspec1,:)) 1.2*max(Output(Viewspec1,:))]);
ax = gca; ax.YAxis.Exponent = 0;

subplot(2,2,4)
plot(EN,Output(Viewspec2,:),'r-');
title(['Spectrum ' num2str(Viewspec2) ':, ChiSq ' num2str(RChiSq(Viewspec2)) ', Median ChiSq ' num2str(MRChiSq)]);
xlim([EN(Windows(1)-PlotBuffer) EN(Windows(4)+PlotBuffer)]);
ylim([1.2*min(Output(Viewspec2,:)) 1.2*max(Output(Viewspec2,:))]);
ax = gca; ax.YAxis.Exponent = 0;

num2str(MaxInd)

% figure('Name','2','units','normalized','outerposition',[0 0 1 1]);
% set(gcf, 'color', [1 1 1]);
% 
% subplot(2,2,1)
% FitSpec = (Input(MaxInd,:));
% fitwin1=FitSpec(Windows(1):Windows(2));
% fitwin2=FitSpec(Windows(3):Windows(4));
% fitwin = cat(2,fitwin1,fitwin2); %Concatenates the values from the two windows. 
% fitEN = cat(2,EN1,EN2);
% %Fit Gaussian plus Lorentzian to the peaks.
% %Fit Gaussian plus Lorentzian to the peaks.
% func= @(C,x)PearsonVII(x,C(1),C(2),C(3),C(4))+C(5);
% C1 = [max(FitSpec),          E0,           FWHM/2, Powguess, median(fitwin2)]; % Initial guesses for the fit parameters
% lb = [max(FitSpec)/2,  E0-E0Var, (FWHM-FWHMVar)/2,   Powmin, min(FitSpec)]; % Lower bounds for the fit parameters
% ub = [max(FitSpec)*40, E0+E0Var, (FWHM+FWHMVar)/2,   Powmax, max(fitwin2)]; % Upper bounds for the fit parameters 
% FR = linspace(fitEN(1),fitEN(end));
% ft = lsqcurvefit(func,C1,fitEN,fitwin, lb, ub, options); % Fit Curve to data.
% ft(2) 
% ft(3) 
% ft(4) 
% ft(5)
% hold on;
% plot(EN,Input(MaxInd,:),'ro', FR,func(ft,FR),'b-');
% title(['Worst Fit Spectrum: ' num2str(MaxInd) ', ChiSq ' num2str(MaxRChiSq)]);
% xlim([EN(Windows(1)-PlotBuffer) EN(Windows(4)+PlotBuffer)]);
% ylim([min(fitwin) max(fitwin)]);
% ax = gca; ax.YAxis.Exponent = 0;
% x1 = EN(Windows(1));
% x2 = EN(Windows(2));
% y2 = [min(Input(MaxInd,:)) max(Input(MaxInd,:))];
% plot([x1 x1],y2, 'b-')
% plot([x2 x2],y2, 'b-')    
% x3 = EN(Windows(3));
% x4 = EN(Windows(4));
% y4 = [min(Input(MaxInd,:)) max(Input(MaxInd,:))];
% plot([x3 x3],y4, 'g-')
% plot([x4 x4],y4, 'g-') 
% box on;
% hold off;
% 
% subplot(2,2,2)
% FitSpec = (Input(MinInd,:));
% fitwin1=FitSpec(Windows(1):Windows(2));
% fitwin2=FitSpec(Windows(3):Windows(4));
% fitwin = cat(2,fitwin1,fitwin2); %Concatenates the values from the two windows. 
% fitEN = cat(2,EN1,EN2);
% %Fit Gaussian plus Lorentzian to the peaks.
% func= @(C,x)PearsonVII(x,C(1),C(2),C(3),C(4))+C(5);
% C1 = [max(FitSpec),          E0,           FWHM/2, Powguess, median(fitwin2)]; % Initial guesses for the fit parameters
% lb = [max(FitSpec)/2,  E0-E0Var, (FWHM-FWHMVar)/2,   Powmin, min(FitSpec)]; % Lower bounds for the fit parameters
% ub = [max(FitSpec)*40, E0+E0Var, (FWHM+FWHMVar)/2,   Powmax, max(fitwin2)]; % Upper bounds for the fit parameters 
%  FR = linspace(fitEN(1),fitEN(end));
% ft = lsqcurvefit(func,C1,fitEN,fitwin, lb, ub, options); % Fit Curve to data.
% hold on;
% plot(EN,Input(MinInd,:),'ro', FR,func(ft,FR),'b-');
% title(['Best Fit Spectrum: ' num2str(MinInd) ', ChiSq ' num2str(MinRChiSq)]);
% xlim([EN(Windows(1)-PlotBuffer) EN(Windows(4)+PlotBuffer)]);
% ylim([min(fitwin) max(fitwin)]);
% ax = gca; ax.YAxis.Exponent = 0;
% x1 = EN(Windows(1));
% x2 = EN(Windows(2));
% y2 = [min(Input(MinInd,:)) max(Input(MinInd,:))];
% plot([x1 x1],y2, 'b-')
% plot([x2 x2],y2, 'b-')    
% x3 = EN(Windows(3));
% x4 = EN(Windows(4));
% y4 = [min(Input(MinInd,:)) max(Input(MinInd,:))];
% plot([x3 x3],y4, 'g-')
% plot([x4 x4],y4, 'g-') 
% box on;
% hold off;
% 
% 
% subplot(2,2,3)
% plot(EN,Output(MaxInd,:),'r-');
% title(['Worst Fit Spectrum: ' num2str(MaxInd) ', ChiSq ' num2str(MaxRChiSq)]);
% xlim([EN(Windows(1)-PlotBuffer) EN(Windows(4)+PlotBuffer)]);
% subplot(2,2,4)
% plot(EN,Output(MinInd,:),'r-');
% title(['Best Fit Spectrum: ' num2str(MinInd) ', ChiSq ' num2str(MinRChiSq)]);
% xlim([EN(Windows(1)-PlotBuffer) EN(Windows(4)+PlotBuffer)]);
% 
Linescan_SubtractP7=single(Output);
WriteTif(Linescan_SubtractP7,String);
