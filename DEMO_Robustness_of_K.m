%----------------------------------------------------------------------------------------------------
% "If you want to find the secrets of the universe, think in terms of energy, frequency & vibration."
%                                                           —— Nikola Tesla
%----------------------------------------------------------------------------------------------------
%  *[Overview]
%   This DEMO is used to demonstrate the robustness of the SRSA method to the K_guess value. 
%----------------------------------------------------------------------------------------------------
%  *[K_guess @ line 76]
%   By assigning different K_guess values, spectral analyses of the same time series can be       ...
%   obtained (as visualized in Figure 2).
%      *[Caption for Figure 2]
%       subplot (a): Normalized Fourier spectrum (baseline reference).
%       subplot (b): 2D SR spectrum using the correct K_true = 4 (as reference).
%       subplot (c): SR spectra obtained using different K_guess values.As discussed in the       ...
%                    manuscript, a higher-dimensional noise subspace increases robustness to the  ...
%                    choice of K_guess. For this reason, subplot (c) uses the maximum L value to  ...
%                    perform the SR analysis. Although the singular spectrum may lose a clear     ...
%                    perception of the characteristic dimension boundary at higher noise levels,  ...
%                    it still provides a useful reference for choosing an initial K_guess. Note   ...
%                    that the difference in the SR spectra obtained with different initial guesses...
%                    is small and hints at the correct K_true value.
%       subplot (d): Frequency-domain stacking of subplot (b). This corresponds to the "stacked   ...
%                    spectrum" described in the manuscript. For simplicity, an averaged spectrum  ...
%                    is shown here.
%       note       ：The noise level can be adjusted at line 34 (recommended range: 0–4).
%----------------------------------------------------------------------------------------------------
% Author      :   Zhifeng Chen
% E-mail      :   zfchen@whu.edu.cn
% Created     :   2026-01-20  
%----------------------------------------------------------------------------------------------------

% Data simulation (complex sequences to display anterograde and retrograde bilateral spectra)
close all;clc;clear;
WNL             =   4;              % Gaussian white noise level control 
K_true          =   4;
%
DeltaT          =   1;              % sampling interval (s)
N               =   1000;           % total number of samples (1)  
t               =   (0:N-1)*DeltaT; % time (s)
%
f1              =   0.02;           % freq1
f2              =   -0.05;          % freq2
f3              =   0.10;           % freq3
f4              =   -0.45;          % freq4
%
y1              =   exp(1i*(2*pi*f1)*t); % com1
y2              =   exp(1i*(2*pi*f2)*t); % com2
y3              =   exp(1i*(2*pi*f3)*t); % com3      
y4              =   exp(1i*(2*pi*f4)*t); % com4 
y               =   y1+y2+y3+y4;         % noiseless time series
noiese          =   normrnd(0,WNL,1,N)+1i*normrnd(0,WNL,1,N);
yn              =   y+noiese;            % y + Gaussian white noise
% Get singular spectrum
P               =   fix(N/2);                           % row window length
L               =   N+1-P;                              % column window length
H               =   zeros(L,P);                         % initialize Hankle matrix
for i = 0:(L-1)
    H(i+1,:) = yn(i+1:i+P);                             % build Hankle matrix
end 
[~,SUM,~]       =   svd(H);                             % perform SVD decomposition  
nonzero         =   nonzeros(SUM);                      
singularV       =   sort(nonzero,'descend');            % singular value            
singularSN      =   1:length(singularV);                % singular value serial number
% Show singular value spectrum
figure(1)
plot(singularSN,singularV, 'd');
% Fourier transform to obtain the Fourier spectrum and normalized for comparison under the same standard.
f0                          =   (0:3*N-1)/DeltaT/N/3;
FFTf                        =   f0-f0(fix(3*N/2)+1);                         % bilateral spectrum
FFTp                        =   abs(fftshift(fft((yn.').*hann(N),3*N)))*2/N; % Fourier transform
NormFFTp                    =   (FFTp-min(FFTp))./(max(FFTp)-min(FFTp));     % Fourier spectrum normalization 
% perform SRSA to obtain SR spectrum
R                           =   1;                   % super-resolution ratio
f0                          =   (0:R*N-1)/DeltaT/N/R;% candidate frequencies for interference-response test  
SRSf                        =   f0-f0(fix(R*N/2)+1);              % bilateral spectrum
K_guess                     =   2; % specify K value, for example, such as pick between [1, 20], even bigger
[SRSf,LVec,SRSp]            =   DoSRSA(yn,DeltaT,SRSf,10,K_guess);% perform SRSA 
[~   ,LVet,SRSt]            =   DoSRSA(yn,DeltaT,SRSf,10,K_true); % perform SRSA (K_true = 4)
%  
figure(2);
tiledlayout(2,2);
% === subplot (a)
nexttile(1)
plot(FFTf,NormFFTp,'k-','LineWidth',1);hold on;
legend('FFT','box','off');
ylabel('Norm Amp');
text(0.02,0.92,'(a)','Units','normalized','FontWeight','bold','FontSize',12);
% === subplot (b)
nexttile(2)
pcolor(SRSf,LVet(:,2)/10,SRSt);
shading interp;colormap('pink');
ylabel('2D SRS');
text(0.02,0.92,'(b)','Units','normalized','FontWeight','bold','FontSize',12,'Color','w');
text(0.2,0.95, sprintf('$K_{ture} = %d$',K_true),'Interpreter','latex', ...
                       'Color','w','Units','normalized','FontSize',12,'FontWeight','bold');
% === subplot (c)               
nexttile(3)
plot(SRSf,SRSp(1,:),'b-','LineWidth',1);
legend('SRSA($L=L_{max}$)','Interpreter','latex','box','off');
ylabel('Norm Amp');
text(0.020,0.92,'(c)','Units','normalized','FontWeight','bold','FontSize',12);
text(0.2,0.95, sprintf('specify $K_{guess} = %d$',K_guess),'Interpreter','latex', ...
                       'Units','normalized','FontSize',12,'FontWeight','bold');
% === subplot (d)
nexttile(4)
plot(SRSf,mean(SRSt(19:39,:),1),'k-');
xlabel('Frequency (Hz)');ylabel('mean spectrum')
text(0.02,0.92,'(d)','Units','normalized','FontWeight','bold','FontSize',12);
text(0.2,0.95, sprintf('spectral stacking'),'Interpreter','latex', ...
                       'Units','normalized','FontSize',12,'FontWeight','bold');
