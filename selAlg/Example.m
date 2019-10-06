clc;clear all
%--------------------------------------------------------------------------------------------------------------
% Recently, I noticed that some papers argued that the BAB-based antenna selection algorithms are sub-optimal. 
% I just worry that they were not using my functions correctly. 
% Therefore, I add this script to give an example how to use my functions correctly. 
% For any questions, please contact me gy519704@outlook.com
% Oct. 06, 2019
%---------------------------------------------------------------------------------------------------------------

N           = 16; 
K           = 4;  
L           = 5;  % suppose an NxK channel matrix, from which we are going to select an optimal LxK submatrix that can achieve maximum channel capcity over all possible submatrices of size LxK. 
rho         = 10; % SNR in linear scale
numTrials   = 10000;

rng(11);
for idxTrial = 1:numTrials
    
    if mod(idxTrial,100)==0
        display(idxTrial);
    end
    
    H = randn(N,K)+1i*randn(N,K); % generate realizations of iid channel matric
    
    % BAB FAS
    % Two things should be noted here when calling BAB based functions 
    % 1. the input channel matrix should be a tall matrix, i.e., # of rows larger than # of columns. If not, make a transpose.
    % 2. the argument SNR rho is the original SNR without normalization to K, since the SNR normalization is already considered inside the function; 
    % 3. the initial bound for BAB-max-capacity algorithm should be -inf; for BAB-min-MSE, the initial bound should be +inf
    selAnt_BAB_FAS      = OptFBB_MaxCap(H.',L,rho,-inf); 
    cap_BAB_FAS         = real(log2(det(eye(K)+rho/K*H(selAnt_BAB_FAS,:)'*H(selAnt_BAB_FAS,:)))); % NOTE, when calculating capcity, the SNR should be normalized

    % Exhaustive search based FAS
    selAnt_ES_FAS       = ESmaxCap(H,L,rho,'fullarray');
    cap_ES_FAS          = real(log2(det(eye(K)+rho/K*H(selAnt_ES_FAS,:)'*H(selAnt_ES_FAS,:))));
    
    assert(isempty(setdiff(selAnt_BAB_FAS, selAnt_ES_FAS))&isempty(setdiff(selAnt_ES_FAS, selAnt_BAB_FAS)),'BAB and ES have different results!!')

end

% toc