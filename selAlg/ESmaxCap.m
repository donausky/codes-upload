% we assume that this is at the Rx side antenna selection
function idxSelectedAnts = ESmaxCap(H,L,rho,flag)

[N,K] = size(H);
if N<K
    error('the number of rows should be larger than the number of columns!');
end
if strcmp(flag,'subarray')
    M           = L;
    numCombs    = M^K;
    idxAll      = permn(1:M,K);
    idxAll      = idxAll + repmat((0:K-1)*M,numCombs,1);
    detVec = zeros(numCombs,1);
    for idxTrial = 1:numCombs
        idxAnts = idxAll(idxTrial,:);
        Hs = H(idxAnts,:);
        detVec(idxTrial) = det(eye(K)+rho/K*Hs'*Hs);%trace(eye(K)/(Hs'*Hs + K/rho*eye(K)));
    end
    [~,idxSelectedComb] = max(detVec);
    idxSelectedAnts = idxAll(idxSelectedComb,:);
elseif strcmp(flag,'fullarray')
    idxAll = nchoosek(1:N,L);
    numCombs = size(idxAll,1);
    detVec = zeros(numCombs,1);
    for idxTrial = 1:numCombs
        idxAnts = idxAll(idxTrial,:);
        Hs = H(idxAnts,:);
        detVec(idxTrial) = det(eye(K)+rho/K*Hs'*Hs);%trace(eye(K)/(Hs'*Hs + K/rho*eye(K)));
    end
    [~,idxSelectedComb] = max(detVec);
    idxSelectedAnts = idxAll(idxSelectedComb,:);
else
    error('undefined array! It should be subarray or fullarray!')
end