function [selAntSet] = GreedyMaxCap(H,Lr,rho,flag)
% this is the greedy antenna selection  to max C AT TX SIDE! If you want to use it at the Rx side, please input H.' instead of H
% Algorithms presented in "Y Gao et al, Massive MIMO Antenna Selection: Switching Architectures, Capacity Bounds and Optimal Antenna Selection Algorithms, IEEE trans. signal processing, 2018"

% Copyright reserved by Yuan Gao, gaoyuan88@gmail.com, San Diego, CA, USA

% input: 
%       H:          channel matrix (a fat matrix). 
%       Lr:         the number of selected antennas
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
G           = rho/Nr*eye(Nr);
ksi         = zeros(Nt,1);
nu          = real(diag(H'*G*H));                   % O(Nu*Nb)
I           = 1:Nt;
if strcmp(flag,'subarray')
    for gg = 1:Lr % loop the groups/blocks
        [~,J]   = max(nu);
        selAntSet(gg)=J;
        groupIdx    = floor((J-1)/M)+1;
        groupElemIdx = (groupIdx-1)*M+1:groupIdx*M;
        I       = setdiff(I,groupElemIdx);
        if gg < Lr
            Delta       = G*H(:,J)/sqrt(1+nu(J));
            G           = G-Delta*Delta';
            ksi(I)      = H(:,I)'*Delta;
            nu(I)       = nu(I) - abs(ksi(I)).^2;
        end
        nu(groupElemIdx)=0;
    end   
elseif strcmp(flag,'fullarray')
    for gg = 1:Lr % loop the groups/blocks
        [~,J]   = max(nu);
        selAntSet(gg)=J;
        I       = setdiff(I,J);
        if gg < Lr
            Delta       = G*H(:,J)/sqrt(1+nu(J));                               % O(Nu^2*Lr)
            G           = G-Delta*Delta';                                       % O(Nu^2*Lr)
            ksi(I)      = H(:,I)'*Delta;                                        % O(Nu*Nb*Lr)
            nu(I)       = nu(I) - abs(ksi(I)).^2;                               % O(Nb*Lr)
        end
        nu(J)=0;
    end   
else
    error('undefined!')
end
% the overall order of complexity is given by O(max{Nu,Nb}*Nu*Lr) = O(NbNuLr)