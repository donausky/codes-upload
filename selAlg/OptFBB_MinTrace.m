function [selAntSet, numNodes] = OptFBB_MinTrace(H,L,rho,initVal)
% this is the full-array branch and bound searching method to search the
% min trace of MSE at the TX side!! If you want to use it at the Rx side, please input H' instead of H
% Algorithms presented in "Y Gao et al, Rotman lens based hybrid analog-digital beamforming in massive MIMO systems: Array architectures beam selection algorithms and experiments, IEEE trans. VT, 2017"

% Copyright reserved by Yuan Gao, gaoyuan88@gmail.com, San Diego, CA, USA

% input: 
%       H:          channel matrix (fat matrix)
%       L:          the number of beams to be selected
%       initVal:    initial Bound, should be large enough! (e.g., positive infinity)
% output:
%       selAntSet:  the selected ant. set
%       numNodes:   number of visited nodes, to show the complexity


[n,m] = size(H);
if m<n
    error('NOTE: this is TX side!! The number of rows should be less than that of columns')
end
% add examine this variable before make it global
global BOUND;
global RETVAL;
global TMPVAL;
global SEARCHFLAG;
global NUMUE
global Nb
global NUMSEL; % number of beams to be selected
global NUMNODES
global Z
global maxIdxEachLevel 

BOUND           = initVal;
NUMSEL          = L;
RETVAL          = zeros(NUMSEL, 1);
TMPVAL          = zeros(NUMSEL, 1);
NUMUE           = n;
Nb              = m;
SEARCHFLAG      = 0;
maxIdxEachLevel = [m-NUMSEL+1:m];
NUMNODES        = 0;

rhoBarInv       = NUMUE/rho;
power_cols      = sum(abs(H).^2,1);                      % O(Nu*Nb)
zeta            = zeros(NUMSEL,1);
for ll = 1:NUMSEL
    zeta(ll)=max(power_cols((1:Nb-NUMSEL+1)+(ll-1)));
end
Z               = zeta./(rhoBarInv+zeta);

G               = eye(NUMUE);
ksi             = zeros(m,1);
A               = H';%H'*G;                               % Note: A and nu are calculated in a recursive form
nu              = power_cols.';%real(diag(H'*G*H));                  
Delta           = nu./(rhoBarInv+nu);                     % O(Nb), the decremental value of Tn 
Dn              = NUMUE;
parentNodeIdx = 0; % start from the root node
optfbb_core(H,1,parentNodeIdx,rhoBarInv,G,A,nu,Dn,ksi,Delta);

if SEARCHFLAG > 0
    selAntSet = RETVAL;%r = RETVAL;
else
    selAntSet = []; %r = 0;
end
numNodes = NUMNODES;

clear BOUND NUMSEL RETVAL TMPVAL NUMUE SEARCHFLAG maxIdxEachLevel;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function optfbb_core(H,layer,parentNodeIdx,rhoBarInv,G,A,nu,Dn,ksi,Delta)

global BOUND;
global RETVAL;
global TMPVAL;
global SEARCHFLAG;
global Nb
global NUMSEL;
global NUMNODES
global Z
global maxIdxEachLevel


if (layer == NUMSEL)
    I               = parentNodeIdx+1:maxIdxEachLevel(NUMSEL);%indices of the child nodes
    DeltaTmp        = Delta(I);
    dunsort         = Dn-DeltaTmp+Z(layer);
    [minVal,minIdx] = min(dunsort);
    if minVal < BOUND
        TMPVAL(NUMSEL)  = I(minIdx);
        RETVAL          = TMPVAL;
        BOUND           = minVal;
        SEARCHFLAG      = SEARCHFLAG + 1;
    end
else
    I           = parentNodeIdx+1:maxIdxEachLevel(layer);     %indices of the child nodes    
    DeltaTmp    = Delta(I);
    dunsort     = Dn-DeltaTmp+Z(layer);
    [~,sortIdx] = sort(dunsort,'ascend');        
    Gtmp        = G;
    DnTmp       = Dn;
    Atmp        = A;
    nuTmp       = nu;
    for ii = sortIdx.'       
        d               = dunsort(ii);
        if (d < BOUND)
            NUMNODES        = NUMNODES+1;
            TMPVAL(layer)   = I(ii);
            Iprime          = I(ii)+1:Nb;
            J               = TMPVAL(layer);
            g               = Gtmp*H(:,J)/sqrt(rhoBarInv+nuTmp(J));
            G               = Gtmp-g*g'; 
            Dn              = DnTmp - DeltaTmp(ii) + Z(layer);            
            ksi(Iprime)     = H(:,Iprime)'*g;
            A(Iprime,:)     = Atmp(Iprime,:) - diag(ksi(Iprime))*repmat(g',length(Iprime),1);
            nu(Iprime)      = nuTmp(Iprime) - abs(ksi(Iprime)).^2;
            Delta(Iprime)   = sum(abs(A(Iprime,:)).^2,2)./(rhoBarInv+nu(Iprime));
            optfbb_core(H,layer+1,I(ii),rhoBarInv,G,A,nu,Dn,ksi,Delta);
        else
            break;
        end
    end
end



