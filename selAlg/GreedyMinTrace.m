function [selAntSet] = GreedyMinTrace(H,Lr,rho,flag)
% this is the greedy antenna selection  to min the trace AT TX SIDE! If you want to use it at the Rx side, please input H' instead of H
% Algorithms presented in "Y Gao et al, Rotman lens based hybrid analog-digital beamforming in massive MIMO systems: Array architectures beam selection algorithms and experiments, IEEE trans. VT, 2017"

% Copyright reserved by Yuan Gao, gaoyuan88@gmail.com, San Diego, CA, USA

% input: 
%       H:          channel matrix (a fat matrix) 
%       Lr:         the number of selected ants
%       rho:        SNR: totalPower/N0 !!
%       flag:       'fullarray','subarray'
% output:
%       selAntSet:  the selected ant. set
% rho = 10;
[Nr,Nt]     = size(H);    
M           = Nt/Lr; % number of elements in a subarray
if Nr>Nt
    error('Channel matrix should be a fat matrix!')
end

selAntSet   = zeros(Lr,1);
% alpha       = sum(abs(H).^2,1).';
G           = rho/Nr*eye(Nr);
ksi         = zeros(Nt,1);
A           = H'*G;                                 % O(Nu^2*Nb)
nu          = real(diag(H'*G*H));                   % O(Nu*Nb)
alpha       = sum(abs(A).^2,2)./(1+nu);             % O(Nu*Nb)                
I           = 1:Nt;
if strcmp(flag,'subarray')
    for gg = 1:Lr % loop the groups/blocks
        [~,J]   = max(alpha);
        selAntSet(gg)=J;
        groupIdx    = floor((J-1)/M)+1;
        groupElemIdx = (groupIdx-1)*M+1:groupIdx*M;
        I       = setdiff(I,groupElemIdx);
        if gg < Lr
            Delta       = G*H(:,J)/sqrt(1+nu(J));
            G           = G-Delta*Delta';
            ksi(I)      = H(:,I)'*Delta;
            A(I,:)      = A(I,:) - diag(ksi(I))*repmat(Delta',length(I),1);
            nu(I)       = nu(I) - abs(ksi(I)).^2;
            alpha(I)    = sum(abs(A(I,:)).^2,2)./(1+nu(I));        
        end
        alpha(groupElemIdx)=0;
    end   
elseif strcmp(flag,'fullarray')
    for gg = 1:Lr % loop the groups/blocks
        [~,J]   = max(alpha);
        selAntSet(gg)=J;
        I       = setdiff(I,J);
        if gg < Lr
            Delta       = G*H(:,J)/sqrt(1+nu(J));                               % O(Nu^2*Lr)
            G           = G-Delta*Delta';                                       % O(Nu^2*Lr)
            ksi(I)      = H(:,I)'*Delta;                                        % O(Nu*Nb*Lr)
            A(I,:)      = A(I,:) - diag(ksi(I))*repmat(Delta',length(I),1);     % O(Nu*Nb*Lr)
            nu(I)       = nu(I) - abs(ksi(I)).^2;                               % O(Nb*Lr)
            alpha(I)    = sum(abs(A(I,:)).^2,2)./(1+nu(I));                     % O(Nu*Nb*Lr)
        end
        alpha(J)=0;
    end   
else
    error('undefined!')
end
% the overall order of complexity is given by O(max{Nu,Nb}*Nu*Lr) = O(NbNuLr)