function [Bins,Q_LVEC,IM] = DoSRSA(CObsY,DeltaT,bins,TPG,varargin)
% DoSRSA Compact SRSA implementation with explicit Q_L/L notation. The  ...
%   code provides two computational modes: full singular value          ...
%   decomposition (full SVD) and randomized singular value decomposition...
%   (rSVD). Full SVD produces deterministic, high-accuracy decomposition...
%   results at a relatively high computational cost, whereas rSVD       ...
%   improves computational efficiency by approximating the dominant     ...
%   singular subspace through randomized projection, with a potential   ...
%   loss of numerical accuracy.
%
% Users may further optimize the computational workflow and             ...
%   parallelization strategy according to the available hardware. When  ...
%   multicore CPUs, GPUs, or multiple computing nodes are available,    ...
%   parallel processing and hardware acceleration can substantially     ...
%   reduce the computational time. The present code prioritizes         ...
%   algorithmic completeness, readability, and reproducibility, and     ...
%   therefore provides only a minimal reference implementation without  ...
%   hardware-specific optimization.
%
%
%   [Bins,Q_LVEC,IM] = DoSRSA(CObsY,DeltaT,bins,TPG)
%   [Bins,Q_LVEC,IM] = DoSRSA(...,'Q_L',Q_L)
%   [Bins,Q_LVEC,IM] = DoSRSA(...,'AlgorithmMode',mode)
%
%   AlgorithmMode:
%       'svd'  : explicit Hankel matrix and complete SVD (default)
%       'rsvd' : matrix-free randomized SVD with complete-SVD fallback
%
%   Notation:
%       Q_L : number of Hankel rows
%       L   : number of Hankel columns/window length, L = N-Q_L+1
%
%   The function automatically samples at most 50 integer Q_L values 
%
%   Q_LVEC columns are [row_index, Q_L, L].

%----------------------------------------------------------------------------------------------------
% Author       :   Zhifeng Chen
% E-mail       :   zfchen@whu.edu.cn
% Created      :   2026-01-20
% Last updated :   2026-06-14
% ==========================================================================
AutoMaxNumQ_L       = 50;
RSVDOversampling    = 10;
RSVDPowerIter       = 1;
RSVDSeed            = 5489;
FreqBlock           = 2048;

CObsY               = CObsY(:);
bins                = bins(:).';

N                   = numel(CObsY);
len                 = numel(bins);

if N < 3
    error('DoSRSA:InputTooShort','CObsY must contain at least 3 samples.');
end

if len < 1
    error('DoSRSA:EmptyBins','bins must contain at least one frequency.');
end

if any(~isfinite(CObsY)) || any(~isfinite(bins))
    error('DoSRSA:InvalidInput','CObsY and bins must contain only finite values.');
end

if ~(isscalar(DeltaT) && isfinite(DeltaT) && DeltaT > 0)
    error('DoSRSA:InvalidDeltaT','DeltaT must be a positive finite scalar.');
end

if ~(isscalar(TPG) && isfinite(TPG) && TPG == round(TPG) && TPG > 0)
    error('DoSRSA:InvalidTPG','TPG must be a positive integer scalar.');
end

if TPG >= N
    error('DoSRSA:TPGTooLarge','TPG must be smaller than N.');
end

p = inputParser;
p.FunctionName = mfilename;
p.PartialMatching = false;
addParameter(p,'Q_L',NaN);
addParameter(p,'AlgorithmMode','svd');
parse(p,varargin{:});

Q_LInput = p.Results.Q_L;
AlgorithmMode = validatestring(lower(char(p.Results.AlgorithmMode)), ...
    {'svd','rsvd'},mfilename,'AlgorithmMode');
isDefaultQL = isnumeric(Q_LInput) && isscalar(Q_LInput) && isnan(Q_LInput);

Q_LStart = 1;
Q_LEnd = N - TPG;

if isDefaultQL
    rangeCount = Q_LEnd - Q_LStart + 1;

    if rangeCount <= AutoMaxNumQ_L
        Q_LV = Q_LStart:Q_LEnd;
    else
        Q_LV = unique(round(linspace( ...
            Q_LStart,Q_LEnd,AutoMaxNumQ_L)),'stable');
    end
else
    if ~isnumeric(Q_LInput) || isempty(Q_LInput)
        error('DoSRSA:InvalidQL','Q_L must be a non-empty numeric vector.');
    end

    Q_LRequested = Q_LInput(:).';

    if any(~isfinite(Q_LRequested)) || ...
            any(Q_LRequested ~= round(Q_LRequested))
        error('DoSRSA:InvalidQL','All Q_L values must be finite integers.');
    end

    Q_LV = unique(round(Q_LRequested),'stable');
    validMask = (Q_LV >= Q_LStart) & (Q_LV <= Q_LEnd);
    numRemoved = sum(~validMask);

    if numRemoved > 0
        warning('DoSRSA:InvalidQLRemoved', ...
                '%d invalid Q_L values are removed. Need 1 <= Q_L <= N-TPG.',...
                 numRemoved);
    end

    Q_LV = Q_LV(validMask);

    if isempty(Q_LV)
        error('DoSRSA:NoValidQL','No valid Q_L remains. Need 1 <= Q_L <= N-TPG.');
    end
end

lenQ_L = numel(Q_LV);
LVec = N - Q_LV + 1;
IM = zeros(lenQ_L,len);
omega = 2*pi*DeltaT*bins;

parfor k = 1:lenQ_L
    Q_L = Q_LV(k);
    L = N - Q_L + 1;

    if strcmp(AlgorithmMode,'svd') || TPG >= min(Q_L,L)
        IM(k,:) = localCompleteSVDRow( ...
            CObsY,Q_L,L,TPG,bins,DeltaT);
    else
        Ps = localRSVDSignalSubspace(CObsY,Q_L,L,TPG, ...
            RSVDOversampling,RSVDPowerIter,RSVDSeed);
        IM(k,:) = localProjectionComplement( ...
            Ps,L,omega,len,FreqBlock);
    end
end

rowMin = min(IM,[],2);
rowMax = max(IM,[],2);
rowRange = rowMax - rowMin;
rowRange(rowRange == 0) = 1;
IM = (IM-rowMin)./rowRange;

if ~isreal(CObsY)
    IM = fliplr(IM);
    Bins = -fliplr(bins);
else
    Bins = bins;
end

Q_LVEC = [(1:numel(Q_LV(:))).',Q_LV(:),LVec(:)];

end


%% ========================================================================================
function IMrow = localCompleteSVDRow(CObsY,Q_L,L,TPG,bins,DeltaT)

HM = hankel(CObsY(1:Q_L),CObsY(Q_L:end));
V = transpose(exp(1i*2*pi*(bins')*DeltaT*(0:L-1)));
[~,~,P] = svd(HM);
Pn = P(:,TPG+1:end);

IMrow = abs(dot(V,V)./(dot((V'*Pn).',(V'*Pn).')));

end

%% ========================================================================================
function Ps = localRSVDSignalSubspace(CObsY,Q_L,L,TPG,oversampling,powerIter,baseSeed)

ops = localHankelOperator(CObsY,Q_L,L);
r = min(TPG+oversampling,min(Q_L,L));

seed = mod(baseSeed + 1009*double(Q_L) + 9173*double(L) + 37*double(TPG), ...
    2^32-1);
stream = RandStream('mt19937ar','Seed',seed);

if isreal(CObsY)
    Omega = randn(stream,L,r);
else
    Omega = complex(randn(stream,L,r),randn(stream,L,r));
end

Y = localHankelMultiply(ops,Omega);
[Qbasis,~] = qr(Y,0);

for ii = 1:powerIter
    Z = localHankelTransposeMultiply(ops,Qbasis);
    [Qz,~] = qr(Z,0);
    Y = localHankelMultiply(ops,Qz);
    [Qbasis,~] = qr(Y,0);
end

B = localHankelTransposeMultiply(ops,Qbasis)';
[~,~,V] = svd(B,'econ');
Ps = V(:,1:TPG);

end

%% ========================================================================================
function ops = localHankelOperator(CObsY,Q_L,L)

N = numel(CObsY);
nfft = 2^nextpow2(max(1,2*N-1));

ops.Q_L = Q_L;
ops.L = L;
ops.IsReal = isreal(CObsY);
ops.NFFT = nfft;
ops.FY = fft(CObsY,nfft);
ops.FCY = fft(conj(CObsY),nfft);

end

%% ========================================================================================
function y = localHankelMultiply(ops,x)

FX = fft(flipud(x),ops.NFFT,1);
C = ifft(ops.FY.*FX,[],1);
y = C(ops.L:(ops.L+ops.Q_L-1),:);

if ops.IsReal && isreal(x)
    y = real(y);
end

end

%% ========================================================================================
function y = localHankelTransposeMultiply(ops,x)

FX = fft(flipud(x),ops.NFFT,1);
C = ifft(ops.FCY.*FX,[],1);
y = C(ops.Q_L:(ops.Q_L+ops.L-1),:);

if ops.IsReal && isreal(x)
    y = real(y);
end

end

%% ========================================================================================
function IMrow = localProjectionComplement(Ps,L,omega,lenFreq,freqBlock)

IMrow = zeros(1,lenFreq);
t = (0:L-1).';
BH = Ps';
Lval = double(L);

for i1 = 1:freqBlock:lenFreq
    i2 = min(i1+freqBlock-1,lenFreq);
    idx = i1:i2;

    Vb = exp(1i*t*omega(idx));
    energy = sum(abs(BH*Vb).^2,1);
    IMrow(idx) = abs(Lval./real(Lval-energy));
end

end
