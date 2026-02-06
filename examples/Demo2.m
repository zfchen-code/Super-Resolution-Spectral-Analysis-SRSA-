%----------------------------------------------------------------------------------------------------
% "If you want to find the secrets of the universe, think in terms of energy, frequency & vibration."
%                                                           —— Nikola Tesla
%----------------------------------------------------------------------------------------------------
%  *[Overview]
%   This script is used to demonstrate the robustness of the SRSA method to the K_used value. 
%----------------------------------------------------------------------------------------------------
%   By assigning different K_used values, spectral analyses of the same time series can be        ...
%   obtained (as visualized in Figure 2).
%      *[Caption for Figure 2]
%       subplot (a): SR spectrum specify K_used=2.
%       subplot (b): fourier spectrum.
%       subplot (c): SR spectrum specify K_used=7.
%       subplot (d): 2D SR spectrum using the K_true = 4 (as reference).
%       subplot (e): SR spectrum specify K_used=10.
%----------------------------------------------------------------------------------------------------
% Author      :   Zhifeng Chen
% E-mail      :   zfchen@whu.edu.cn
% Created     :   2026-01-20  
%----------------------------------------------------------------------------------------------------
% Data simulation (complex sequences to display anterograde and retrograde bilateral spectra)
clc;clear;
%
DeltaT          =   1;              % sampling interval (s)
N               =   1000;           % total number of samples (1)  
t               =   (0:N-1)*DeltaT; % time (s)
%                  AMP              Freq
yn              =   1*exp(1i*(2*pi* 0.02)*t)+...
                    1*exp(1i*(2*pi*-0.05)*t)+... 
                    1*exp(1i*(2*pi* 0.10)*t)+...      
                    1*exp(1i*(2*pi*-0.45)*t)+normrnd(0,4,1,N)+1i*normrnd(0,4,1,N);                         
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
% Figure 1 
fig4=figure(4);
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
K_used1                     =   2; % specify K value, for example, such as pick between [1, 20], even bigger
[SRSf,LVec1,SRSp1]          =   DoSRSA(yn,DeltaT,SRSf,10,K_used1);% perform SRSA 
K_used2                     =   7; % specify K value, for example, such as pick between [1, 20], even bigger
[SRSf,LVec2,SRSp2]          =   DoSRSA(yn,DeltaT,SRSf,10,K_used2);% perform SRSA 
K_used3                     =   10;% specify K value, for example, such as pick between [1, 20], even bigger
[SRSf,LVec3,SRSp3]          =   DoSRSA(yn,DeltaT,SRSf,10,K_used3);% perform SRSA 
K_true                      =   4;
[~   ,LVet,SRSt]            =   DoSRSA(yn,DeltaT,SRSf,10,K_true); % perform SRSA (K_true = 4)
% Figure 2 
figure(2);tiledlayout(6,2);
nexttile(1,[2,1])
plot(SRSf,SRSp1(1,:),'k-','LineWidth',1);
ylim([0,1.2]);
legend('SRSA($L=L_{max}$)','Interpreter','latex','box','off');
ylabel('Norm Amp');
text(0.020,0.92,'(a)','Units','normalized','FontWeight','bold','FontSize',12);
text(0.2,0.92, sprintf('specify $K_{used} = %d$',K_used1),'Interpreter','latex', ...
                       'Units','normalized','FontSize',12,'FontWeight','bold');
nexttile(2,[3,1])
plot(FFTf,NormFFTp,'k-','LineWidth',1);
ylim([0,1.2]);
ylabel('Norm Amp');
text(0.02,0.95,'(b)','Units','normalized','FontWeight','bold','FontSize',12);
legend('Fourier Spectrum','box','off');
nexttile(8,[3,1])
pcolor(SRSf,LVet(:,2)/10,SRSt);
shading interp;colormap('pink');
ylabel('2D SRS');
yticklabels([]);
text(0.02,0.95,'(d)','Units','normalized','FontWeight','bold','FontSize',12,'Color','w');
text(0.2,0.92, sprintf('$K_{ture} = %d$',K_true),'Interpreter','latex', ...
                       'Color','w','Units','normalized','FontSize',12,'FontWeight','bold');     
legend('2D SRS','box','off','textcolor','w');
nexttile(5,[2,1])
plot(SRSf,SRSp2(1,:),'k-','LineWidth',1);
ylim([0,1.2]);
legend('SRSA($L=L_{max}$)','Interpreter','latex','box','off');
ylabel('Norm Amp');
text(0.020,0.92,'(c)','Units','normalized','FontWeight','bold','FontSize',12);
text(0.2,0.95, sprintf('specify $K_{used} = %d$',K_used2),'Interpreter','latex', ...
                       'Units','normalized','FontSize',12,'FontWeight','bold');
nexttile(9,[2,1])
plot(SRSf,SRSp3(1,:),'k-','LineWidth',1);
ylim([0,1.2]);
legend('SRSA($L=L_{max}$)','Interpreter','latex','box','off');
ylabel('Norm Amp');
text(0.020,0.92,'(e)','Units','normalized','FontWeight','bold','FontSize',12);
text(0.2,0.92, sprintf('specify $K_{used} = %d$',K_used3),'Interpreter','latex', ...
                       'Units','normalized','FontSize',12,'FontWeight','bold');
