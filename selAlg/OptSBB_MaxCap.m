function [selAntSet, numIters] = OptSBB_MaxCap(H,M,rho,initVal)
% this is the sub-array branch and bound searching method to search the
% max. channel capacity at the TX side!! If you want to use it at the Rx side, please input H' instead of H
% Algorithms presented in "Y Gao et al, Massive MIMO Antenna Selection: Switching Architectures, Capacity Bounds and Optimal Antenna Selection Algorithms, IEEE trans. signal processing, 2018"

% input: 
%       H:          channel matrix
%       M:          subblocksize
%       initVal:    initial Bound, should be small enough! (e.g., negative infinity)
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

rhoBar          = rho/NUMUE;
power_cols      = sum(abs(H).^2,1);   
zeta            = max(reshape(power_cols,M,[]),[],1).';
Z               = log2(1+rhoBar*zeta);

G               = eye(NUMUE);
ksi             = zeros(m,1);
nu              = power_cols.';%real(diag(H'*G*H));                  
Delta           = log2(1+rhoBar*nu);%nu./(rhoBarInv+nu);                     % O(Nb), the decremental value of Tn 
Dn              = 0;
optsbb_core(H,1,rhoBar,G,nu,Dn,ksi,Delta);

if SEARCHFLAG > 0
    selAntSet = RETVAL;
else
    selAntSet = []; 
end
numIters = NUMITERS;

clear BOUND NUMSB RETVAL TMPVAL NUMUE SUBBLOCKSIZE SEARCHFLAG NUMITERS;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optsbb_core(H,layer,rhoBar,G,nu,Dn,ksi,Delta)

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
    DeltaTmp            = Delta(I);%log2(1+rhoBar*Delta(I));
    dunsort             = Dn+DeltaTmp-Z(layer);
    [maxVal,maxIdx]     = max(dunsort);
    if maxVal > BOUND
        TMPVAL(NUMSB)   = maxIdx + (NUMSB-1)*SUBBLOCKSIZE;
        RETVAL          = TMPVAL;
        BOUND           = maxVal;
        SEARCHFLAG      = SEARCHFLAG + 1;
    end
else
    I                   = (1:SUBBLOCKSIZE)+(layer-1)*SUBBLOCKSIZE;    
    DeltaTmp            = Delta(I);
    dunsort             = Dn+DeltaTmp-Z(layer);
    [~,sortIdx]         = sort(dunsort,'descend');        
    Gtmp                = G; % backup these variables are very important due to the FOR LOOP
    DnTmp               = Dn;
    nuTmp               = nu;
    for ii = sortIdx.'        
        d               = dunsort(ii);
        if (d > BOUND)
            NUMITERS        = NUMITERS+1;
            TMPVAL(layer)   = ii+(layer-1)*SUBBLOCKSIZE;
            Iprime          = layer*SUBBLOCKSIZE+1:NUMSB*SUBBLOCKSIZE; % the residual beam indecies need to be updated
            J               = TMPVAL(layer);
            g               = Gtmp*H(:,J)/sqrt(1/rhoBar+nuTmp(J));
            G               = Gtmp-g*g'; 
            Dn              = DnTmp + DeltaTmp(ii) - Z(layer);            
            ksi(Iprime)     = H(:,Iprime)'*g;
            nu(Iprime)      = nuTmp(Iprime) - abs(ksi(Iprime)).^2;
            Delta(Iprime)   = log2(1+rhoBar*nu(Iprime));%sum(abs(A(Iprime,:)).^2,2)./(rhoBarInv+nu(Iprime));
            optsbb_core(H,layer+1,rhoBar,G,nu,Dn,ksi,Delta);
        else
            break;
        end
    end
end



