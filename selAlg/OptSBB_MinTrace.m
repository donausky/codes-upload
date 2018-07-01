function [selAntSet, numIters] = OptSBB_MinTrace(H,M,rho,initVal)
% this is the sub-array branch and bound searching method to search the
% min trace of MSE at the TX side!! If you want to use it at the Rx side, please input H.' instead of H
% Algorithms presented in "Y Gao et al, Rotman lens based hybrid analog-digital beamforming in massive MIMO systems: Array architectures beam selection algorithms and experiments, IEEE trans. VT, 2017"

% Copyright reserved by Yuan Gao, gaoyuan88@gmail.com, San Diego, CA, USA

% input: 
%       H:          channel matrix (a fat matrix)
%       M:          subblocksize
%       initVal:    initial Bound, should be large enough! (e.g., positive infinity)
% output:
%       selAntSet:  the selected ant. set
%       numIters:   number of visited nodes, to show the complexity

[n,m] = size(H);
if m<n
    error('NOTE: this is TX side!! The number of rows should be less than that of columns')
end
% add examine this variable before make it global
global BOUND;
global RETVAL;
global TMPVAL;
global SUBBLOCKSIZE;
global SEARCHFLAG;
global NUMUE
global NUMSB
global NUMITERS
global Z

BOUND           = initVal;
NUMSB           = m/M;
RETVAL          = zeros(NUMSB, 1);
TMPVAL          = zeros(NUMSB, 1);
NUMUE           = n;
SUBBLOCKSIZE    = M;
SEARCHFLAG      = 0;
NUMITERS        = 0;

rhoBarInv       = NUMUE/rho;
power_cols      = sum(abs(H).^2,1);   
zeta            = max(reshape(power_cols,M,[]),[],1).';
Z               = zeta./(rhoBarInv+zeta);%log2(1+rho/NUMUE*zeta);

G               = eye(NUMUE);
ksi             = zeros(m,1);
A               = H';%H'*G;                               % Note: A and nu are calculated in a recursive form
nu              = power_cols.';%real(diag(H'*G*H));                  
Delta           = nu./(rhoBarInv+nu);                     % O(Nb), the decremental value of Tn 
Dn              = NUMUE;
optsbb_core(H,1,rhoBarInv,G,A,nu,Dn,ksi,Delta);

if SEARCHFLAG > 0
    selAntSet = RETVAL;
else
    selAntSet = []; 
end
numIters = NUMITERS;

clear BOUND NUMSB RETVAL TMPVAL NUMUE SUBBLOCKSIZE SEARCHFLAG NUMITERS;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optsbb_core(H,layer,rhoBarInv,G,A,nu,Dn,ksi,Delta)

global BOUND;
global RETVAL;
global TMPVAL;
global SUBBLOCKSIZE;
global SEARCHFLAG;
global NUMUE
global NUMSB
global NUMITERS
global Z

if (layer == NUMSB)
    I                   = (1:SUBBLOCKSIZE)+(layer-1)*SUBBLOCKSIZE;
    DeltaTmp            = Delta(I);%log2(1+rho/NUMUE*Delta(I));
    dunsort             = Dn-DeltaTmp+Z(layer);
    [minVal,minIdx]     = min(dunsort);
    if minVal < BOUND
        TMPVAL(NUMSB)   = minIdx + (NUMSB-1)*SUBBLOCKSIZE;
        RETVAL          = TMPVAL;
        BOUND           = minVal;
        SEARCHFLAG      = SEARCHFLAG + 1;
    end
else
    I                   = (1:SUBBLOCKSIZE)+(layer-1)*SUBBLOCKSIZE;    
    DeltaTmp            = Delta(I);
    dunsort             = Dn-DeltaTmp+Z(layer);
    [~,sortIdx]         = sort(dunsort,'ascend');        
    Gtmp                = G; % backup these variables are very important due to the FOR LOOP
    DnTmp               = Dn;
    Atmp                = A;
    nuTmp               = nu;
    for ii = sortIdx.'        
        d               = dunsort(ii);
        if (d < BOUND)
            NUMITERS        = NUMITERS+1;
            TMPVAL(layer)   = ii+(layer-1)*SUBBLOCKSIZE;
            Iprime          = layer*SUBBLOCKSIZE+1:NUMSB*SUBBLOCKSIZE; % the residual beam indecies need to be updated
            J               = TMPVAL(layer);
            g               = Gtmp*H(:,J)/sqrt(rhoBarInv+nuTmp(J));
            G               = Gtmp-g*g'; 
            Dn              = DnTmp - DeltaTmp(ii) + Z(layer);            
            ksi(Iprime)     = H(:,Iprime)'*g;
            A(Iprime,:)     = Atmp(Iprime,:) - diag(ksi(Iprime))*repmat(g',length(Iprime),1);
            nu(Iprime)      = nuTmp(Iprime) - abs(ksi(Iprime)).^2;
            Delta(Iprime)   = sum(abs(A(Iprime,:)).^2,2)./(rhoBarInv+nu(Iprime));
            optsbb_core(H,layer+1,rhoBarInv,G,A,nu,Dn,ksi,Delta);
        else
            break;
        end
    end
end



