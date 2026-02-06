% This script demonstrates the analysis of real data. The M1 station record from  ...
%   the 2011 event (The m = 0 0S2 component is generally considered difficult to  ...
%   excite to observable levels) is used, and the 0S2 mode is selected as the     ...
%   analysis target to reproduce Fig. 5(h) of the manuscript. Because SRSA is     ...
%   computationally intensive, only a single time series is provided for          ...
%   demonstration. After execution, please allow a short time for the results to  ...
%   be generated.
% *The raw data were preprocessed and corrected for modeled signals to obtain time...
%   series y_1; The time series y_1 was further band-pass filtered (0.2–0.7 mHz)  ...
%   and downsampled to obtain y_2.
%   *[Caption for Figure 3]
%       subplot (a): time series y_1 & time series y_2.
%       subplot (b): Power spectra of y_1 & y_2.
%       subplot (c): Zoom of the power spectrum of y_1 around the 0S2 band.
%       subplot (d): Zoom of the power spectrum of y_2 around the 0S2 band.
%       subplot (e): full singular-value spectrum.
%       subplot (f): Zoom of the elbow region of the singular-value spectrum.
%       subplot (g): 2D SRS obtained with K_used = 50.
%       subplot (h): 2D SRS obtained with K_used = 55.
%       subplot (i): 2D SRS obtained with K_used = 60.
%       subplot (j): mean stacking spectrum of subplot (g).
%       subplot (k): mean stacking spectrum of subplot (h).
%       subplot (l): mean stacking spectrum of subplot (i).
% ----------------------------------------------------------------------------------
clc;clear;
load series_m1.mat
y_1                         =   series_m1{1}(:); % y_1  f_s = 1/60
y_2                         =   series_m1{2}(:); % y_2  f_s = 1/300
%
Win1                        =   hann(numel(y_1),'periodic');
nfft1                       =   2*numel(y_1);
[Power1(:,2),Power1(:,1)]   =   periodogram(y_1,Win1,nfft1,1/60,'power');
seg1                        =   Power1(Power1(:,1)>=0.29e-3 & Power1(:,1)<=0.33e-3,:);
Win2                        =   hann(numel(y_2),'periodic');
nfft2                       =   2*numel(y_2);
[Power2(:,2),Power2(:,1)]   =   periodogram(y_2,Win2,nfft2,1/300,'power');
seg2                        =   Power2(Power2(:,1)>=0.29e-3 & Power2(:,1)<=0.33e-3,:);
%
N1                          =   numel(y_2);
P                           =   fix(N1/2);                 % row window length
L                           =   N1+1-P;                    % column window length
H                           =   zeros(L,P);                % initialize Hankle matrix
for i = 0:(L-1)
    H(i+1,:)                =   y_2(i+1:i+P);              % build Hankle matrix
end 
[~,SUM,~]                   =   svd(H);                    % perform SVD decomposition  
nonzero                     =   nonzeros(SUM);                      
singularV                   =   sort(nonzero,'descend');   % singular value            
singularSN                  =   1:length(singularV);       % singular value serial number
%
Delta_f                     =   1/(300*numel(y_2)*10);     % tenfold super-resolution
bins                        =   0.270e-3:Delta_f:0.350e-3;
[SRSf1,LVec1,SRSp1]         =   DoSRSA(y_2,300,bins,20,50);% perform SRSA K_used = 50
[SRSf2,LVec2,SRSp2]         =   DoSRSA(y_2,300,bins,20,55);% perform SRSA K_used = 55
[SRSf3,LVec3,SRSp3]         =   DoSRSA(y_2,300,bins,20,60);% perform SRSA K_used = 60
%%
figure(3);tiledlayout(6,5)
nexttile(1,[2,2])
plot(1:numel(y_1),y_1,'k-');              hold on
plot(downsample(1:numel(y_1),5),y_2,'r-');hold off
text(0.06,0.92,'(a)','Units','normalized','FontWeight','bold','FontSize',12);
legend('The preprocessed and model-corrected superconducting gravimeter time series, denoted as $y_1$', ...
       'The band-pass-filtered (0.2-0.7 mHz) and downsampled version of $y_1$, denoted as $y_2$', ...
       'Interpreter','latex','box','off','FontSize',10);
ylabel('Amplitude (1)');xlabel('time (min)');
nexttile(11,[2,2])
plot(Power1(:,1)*1000,Power1(:,2),'k-'); hold on;
plot(Power2(:,1)*1000,Power2(:,2),'r--');hold off;
legend('Power spectra of $y_1$', ...
       'Power spectra of $y_2$', ...
       'Interpreter','latex','box','off','FontSize',10);
xlim([0.1 1]);
text(0.06,0.92,'(b)','Units','normalized','FontWeight','bold','FontSize',12);
xlabel('Frequency (mHz)');ylabel('Power Spectra (nm/s^2)^2')
nexttile(21,[2,1])
plot(seg1(:,1)*1000,rescale(seg1(:,2), 0, 1),'k-'); 
xline([0.299965;0.304536;0.309193;0.313843;0.318433],'--','Color',[0.7 0.7 0.7]);
lbl = "$y_1$ Power spectrum in the ${}_0S_2$ band" + newline + "$-\,-\-$ Ref. Ding \& Shen (2013)";
legend(lbl,'Interpreter','latex','FontSize',10,'box','off');
xlim([0.29 0.33]);ylim([0,1.2]);
text(0.06,0.92,'(c)','Units','normalized','FontWeight','bold','FontSize',12);
xlabel('Frequency (mHz)');ylabel('Norm Power Spectra (nm/s^2)^2')
nexttile(22,[2,1])
plot(seg2(:,1)*1000,rescale(seg2(:,2), 0, 1),'r-');
xline([0.299965;0.304536;0.309193;0.313843;0.318433],'--','Color',[0.7 0.7 0.7]);
lbl = "$y_2$ Power spectrum in the ${}_0S_2$ band" + newline + "$-\,-\-$ Ref. Ding \& Shen (2013)";
legend(lbl,'Interpreter','latex','FontSize',10,'box','off');
xlim([0.29 0.33]);ylim([0,1.2]);
text(0.06,0.92,'(d)','Units','normalized','FontWeight','bold','FontSize',12);
xlabel('Frequency (mHz)');ylabel('Norm Power Spectra (nm/s^2)^2')
nexttile(3,[3,1])
plot(singularSN,singularV, 'd','MarkerEdgeColor','k');
legend('Singular-value spectrum of $y_2$','Interpreter','latex','FontSize',10,'box','off')
text(0.06,0.92,'(e)','Units','normalized','FontWeight','bold','FontSize',12);
ylabel('singular value');
nexttile(18,[3,1])
plot(singularSN,singularV, 'd','MarkerEdgeColor','k');hold on
plot(50,singularV(50),'d','MarkerFaceColor','b','MarkerEdgeColor','b','MarkerSize',8);
plot(54,singularV(55),'d','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',8);
plot(60,singularV(60),'d','MarkerFaceColor','g','MarkerEdgeColor','g','MarkerSize',8);hold off
% legend({'Zoom of the elbow region of \n the singular-value spectrum'},'Interpreter','latex','FontSize',10,'box','off')
lbl = "Zoom of the elbow region of the" + newline + "singular-value spectrum subplot (e)";
legend(lbl,'$K_{usd}=50$','$K_{usd}=55$','$K_{usd}=60$','Interpreter','latex','FontSize',10,'box','off');
xlim([30,150]);ylim([0,100]);
text(0.06,0.92,'(f)','Units','normalized','FontWeight','bold','FontSize',12);
xlabel('Series Number');ylabel('singular value');
nexttile(4,[2,1])
pcolor(SRSf1*1000,LVec1(33:117,2),SRSp1(33:117,:));
legend('$K_{usd}=50$','Interpreter','latex','FontSize',10,'box','off','TextColor','b');
xlim([0.29,0.33])
shading interp;colormap('pink');
ylabel('2D SRS');
yticklabels([]);
text(0.06,0.92,'(g)','Units','normalized','FontWeight','bold','FontSize',12,'Color','w');
nexttile(14,[2,1])
pcolor(SRSf2*1000,LVec2(33:117,2),SRSp2(33:117,:));
legend('$K_{usd}=55$','Interpreter','latex','FontSize',10,'box','off','TextColor','r');
xlim([0.29,0.33])
shading interp;colormap('pink');
ylabel('2D SRS');
yticklabels([]);
text(0.06,0.92,'(h)','Units','normalized','FontWeight','bold','FontSize',12,'Color','w');
nexttile(24,[2,1])
pcolor(SRSf3*1000,LVec3(53:117,2),SRSp3(53:117,:));
legend('$K_{usd}=60$','Interpreter','latex','FontSize',10,'box','off','TextColor','g');
xlim([0.29,0.33])
shading interp;colormap('pink');
xlabel('Frequency (mHz)');ylabel('2D SRS');
yticklabels([]);
text(0.06,0.92,'(i)','Units','normalized','FontWeight','bold','FontSize',12,'Color','w');
nexttile(5,[2,1])
plot(SRSf1*1000,mean(SRSp1(33:117,:),1),'b-');
xline([0.299965;0.304536;0.309193;0.313843;0.318433],'--','Color',[0.7 0.7 0.7]);
lbl = "mean stacking spectrum ($K_{usd}=50$)" + newline + "$-\,-\-$ Ref. Ding \& Shen (2013)";
legend(lbl,'Interpreter','latex','FontSize',10,'box','off');
xlim([0.29,0.33])
ylabel('SR Spectrum');
text(0.06,0.92,'(j)','Units','normalized','FontWeight','bold','FontSize',12);
nexttile(15,[2,1])
plot(SRSf2*1000,mean(SRSp2(33:117,:),1),'r-');
xline([0.299965;0.304536;0.309193;0.313843;0.318433],'--','Color',[0.7 0.7 0.7]);
lbl = "mean stacking spectrum ($K_{usd}=55$)" + newline + "$-\,-\-$ Ref. Ding \& Shen (2013)";
legend(lbl,'Interpreter','latex','FontSize',10,'box','off');
xlim([0.29,0.33])
ylabel('SR Spectrum');
text(0.06,0.92,'(k)','Units','normalized','FontWeight','bold','FontSize',12);
nexttile(25,[2,1])
plot(SRSf3*1000,mean(SRSp3(53:117,:),1),'g-');
xline([0.299965;0.304536;0.309193;0.313843;0.318433],'--','Color',[0.7 0.7 0.7]);
lbl = "mean stacking spectrum ($K_{usd}=60$)" + newline + "$-\,-\-$ Ref. Ding \& Shen (2013)";
legend(lbl,'Interpreter','latex','FontSize',10,'box','off');
xlim([0.29,0.33])
xlabel('Frequency (mHz)');ylabel('SR Spectrum');
text(0.06,0.92,'(l)','Units','normalized','FontWeight','bold','FontSize',12);