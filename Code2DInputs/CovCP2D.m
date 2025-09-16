
function K = CovCP2D(theta,X1,X2,SeqX1,SeqX2)
    NStream = [length(SeqX1) length(SeqX2)];
    SS = sum(NStream);
    D = size(X1,2);
    if length(theta) == (D+1)*SS+1
        Hyper.Rho = theta(1:SS);
        Hyper.lambda = theta(SS+(1:D*SS));
    else
        Hyper.Rho = repmat(theta(1:NStream(1)),2,1);
        Hyper.lambda = repmat(theta(NStream(1)+(1:D*NStream(1))),2,1);
    end
    K = CovSub(X1,X2,SeqX1,SeqX2,Hyper,NStream);

    [m,n] = size(K);
    if m == n && sum((X1-X2).^2,"all") == 0
        sigma2 = repmat(theta(end).^2',NStream(1),1);
        sigma2 = repelem(sigma2,SeqX1);
        sigma2 = diag(sigma2);
        K = K+sigma2;
    end

end

function K = CovSub(X1,X2,SeqX1,SeqX2,Hyper,NStream)

rho_ij = Hyper.Rho; % [rho(d1)_i0 rho(d1)_ii rho(d1)_i'0 rho(d1)_i'i';
                    %  rho(d2)_i0 rho(d2)_ii rho(d2)_i'0 rho(d2)_i'i']
                    %  d1/d2--domain index
lambda_ij = (Hyper.lambda).^2; %(exp(Hyper.lambda)).^2; % [l(d1)_i0 l(d1)_ii l(d1)_i'0 l(d1)_i'i';
                          %  l(d2)_i0 l(d2)_ii l(d2)_i'0 l(d2)_i'i']
                          %  d1/d2--domain index

%% Prepare scale parameter

Lm = NStream(1);
Ln = NStream(2);
D = size(X1,2);
RhoProd = rho_ij(1:Lm).*rho_ij(Lm+1:sum(NStream))';
lambda_temp = reshape(lambda_ij,D,[])';
LambdaProd = prod(lambda_temp,2);
LambdaProd = (LambdaProd(1:Lm).*LambdaProd(Lm+1:sum(NStream))').^(1/4);
LambdaSum = zeros(sum(NStream),D);
for i = 1:Ln
    for j = 1:Lm
        LambdaSum(j+(i-1)*Ln,:) = lambda_temp(j,:)+lambda_temp(Lm+i,:);
    end
end
LambdaSumSqrt = sqrt(reshape(prod(LambdaSum,2),Lm,Ln));

RhoProd = repelem(RhoProd,SeqX1,SeqX2);
LambdaProd = repelem(LambdaProd,SeqX1,SeqX2);
LambdaSumSqrt = repelem(LambdaSumSqrt,SeqX1,SeqX2);
Scale = 2^(D/2).*RhoProd.*LambdaProd./LambdaSumSqrt;
%% Prepare the input

X1 = X1./repelem(lambda_temp(1:Lm,:),SeqX1,1);
X2 = X2./repelem(lambda_temp(Lm+1:end,:),SeqX2,1);

D2XX = sq_dist(X1',X2');
K = Scale.*exp(-1/2*D2XX);

end






