function [Bins,LVEC,IM] = DoSRSA(CObsY,DeltaT,bins,LS,TPG,varargin)
%-----------------------------------------------------------------------------------------
%   SRSA：Super-Resolution Spectral Analysis
%   ---------------------------------------
%   [Input]：
%       CObsY           ：time series, real or complex support
%       DeltaT          ：sampling interval
%       bins            ：candidate frequency vector
%       LS              ：window length step
%       TPG             ：initial guess of number of poles
%   [Output]：
%       Bins            ：output frequency vector
%       LVEC            ：window length search vector
%       IM              ：the 2D spectrum output("Interference" Matrix)
%       -----------------------------------
%       Author          : Zhifeng Chen
%       E-mail          : zfchen@whu.edu.cn
%       Created         : 2026-01-20
%       Test platform   : matlab2024b
%-----------------------------------------------------------------------------------------
N           =   numel(CObsY); 
len         =   numel(bins); 
p           =   inputParser;addParameter(p,'L',NaN);parse(p,varargin{:});
if isnan(p.Results.L)
    LVEC    =   1:LS:(N-TPG);
    lenL    =   length(LVEC);
else
    LVEC    =   p.Results.L;
    lenL    =   length(p.Results.L);
end
IM          =   zeros(lenL,len);warning off;
parfor k = 1:lenL % Please develop the best parallel solution based on the local platform
    HM      =   hankel(CObsY(1:LVEC(k)),CObsY(LVEC(k):end));
    V       =   transpose(exp(1i*2*pi*(bins')*DeltaT*(0:(N-LVEC(k)+1)-1))); % \mathbf{V}
    [~,~,P] =   svd(HM);      
    IM(k,:) =   abs(dot(V,V)./(dot((V'*P(:,TPG+1:end)).',(V'*P(:,TPG+1:end)).')));
end
IM          =   (IM-min(IM,[],2))./(max(IM,[],2)-min(IM,[],2));          % normalization
if ~isreal(CObsY);IM = fliplr(IM);Bins = -fliplr(bins);else;Bins = bins;end
LVEC        =   [(1:numel(LVEC'))',LVEC'];
end
