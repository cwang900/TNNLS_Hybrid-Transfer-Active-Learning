
function K = CovLMC2D(theta,X1,X2,SeqX1,SeqX2)
    
    D = size(X1,2);
    NStream = [length(SeqX1) length(SeqX2)];
    NLatent = (length(theta)-1)./(sum(NStream)+D); % Number of latent processes
    Hyper.Rho = theta(1:NLatent*sum(NStream),:)';
    Hyper.lambda = theta(NLatent*sum(NStream)+(1:D*NLatent),:);
    K = CovSub(X1,X2,SeqX1,SeqX2,Hyper,NLatent);

    [m,n] = size(K);
    if m == n && sum((X1-X2).^2,"all") == 0
        sigma2 = theta(end,1).^2';
        sigma2 = repmat(sigma2,sum(SeqX1),1);
        sigma2 = diag(sigma2);
        K = K+sigma2;
    end

end

function K = CovSub(X1,X2,SeqX1,SeqX2,Hyper,NLatent)

D = size(X1,2); % input dimension
NStream = [length(SeqX1) length(SeqX2)];
rho_mm = reshape(Hyper.Rho(:),NLatent,[]);

K = 0;
scale_mm = abs(Hyper.lambda);
scale_mm = reshape(scale_mm(:),D,[]);
for i = 1:NLatent
    
    scale_temp = repmat(scale_mm(:,i),1,NStream(1))';
    X1scale_mm = repelem(scale_temp,SeqX1,1);
    scale_temp = repmat(scale_mm(:,i),1,NStream(2))';
    X2scale_mm = repelem(scale_temp,SeqX2,1);
    X1_mm = X1./X1scale_mm;
    X2_mm = X2./X2scale_mm; 

    D2mm = sq_dist(X1_mm', X2_mm');
    rho_temp = rho_mm(i,:);
    rho_mm_square = rho_temp(1:NStream(1))'.*rho_temp(NStream(1)+1:end);
    rho_mm_square = repelem(rho_mm_square,SeqX1,SeqX2);

    K = K+rho_mm_square.*exp(-1/2*D2mm);
end

end






