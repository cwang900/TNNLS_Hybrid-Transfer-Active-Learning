function [Mu, Sigma, alphat, Ct] = PostInd2D(theta,y,Xf,XStar,alpha,C,XInd,NDomain,CovType,NStream)

theta = theta';

NInd = length(XInd{1});
NT = length(y{end});

%%
if isempty(alpha)
    [KTarget, KSTInd, KTTInd, Gxx_n_n, Inv_KS, A] =...
        CovCalInd2D(theta,y,Xf,XInd,CovType,NDomain,NStream);
    Ainv = eye(NT)/A;

    % \Gamma_sx^{nn}/A
    B = KTTInd'*Ainv;

    % \Gamma_sx^{n,<n}*G_xx^{<n,n}
    Prod = @(K,G) K'*G';
    R = cellfun(Prod,KSTInd(1:NDomain-1),Gxx_n_n(1:NDomain-1),'UniformOutput',false);
    R = cell2mat(R);
    R = reshape(R,NInd,NT,[]);
    R = sum(R,3);

    % G_sx^{n,<n}=\Gamma_sx^{n,<n}*\Sigma_xx^{<n,<n}
    Prod = @(K,G) K'*G;
    Gsx_n_n = cellfun(Prod,KSTInd(1:NDomain-1),Inv_KS,'UniformOutput',false);

    % R/A
    RAinv = R*Ainv;
    B_RAinv = B-RAinv;

    Prod = @(G1,G2,y) (-B_RAinv*G1+G2)*y; 
    if length(y)>1
        MuS = cellfun(Prod,Gxx_n_n,Gsx_n_n,y(1:NDomain-1),'UniformOutput',false);
        MuS = cell2mat(MuS);
        MuS = sum(MuS,2);
    else
        MuS = zeros(NInd,1);
    end
    MuT = B_RAinv*y{end};
    alphat = MuS+MuT;

    RB = R*B';
    Prod = @(G,K) G*K;
    GK = cellfun(Prod,Gsx_n_n,KSTInd(1:NDomain-1),'UniformOutput',false);
    GK = cell2mat(GK);
    GK = reshape(GK,NInd,NInd,[]);
    GK = sum(GK,3);
    Ct = KSTInd{end}-(B*KTTInd-(RB+RB')+GK+RAinv*R')+1e-6*eye(NInd);
else
    [Gxs_nn, A] =...
        CovCalIndOnline2D(theta,y,Xf,XInd,CovType,NDomain,NStream);
    PInv = eye(NInd)/C;
    GA = Gxs_nn'/A;
    Ct = (PInv+GA*Gxs_nn)\eye(NInd)+1e-6*eye(NInd);
    alphat = alpha+C*Gxs_nn'/(A+Gxs_nn*C*Gxs_nn')*(y{end}-Gxs_nn*alpha);%Ct*(GA*y{end}+PInv*alpha);
end

%%

if ~isempty(XStar)
    [Gxs_nn, Qxx_nn, KTStar] =...
        CovCalPred2D(theta,XInd,XStar,NDomain,CovType,NStream);
    Mu = Gxs_nn*alphat;
    Sigma = KTStar-Qxx_nn+Gxs_nn*Ct*Gxs_nn';
else
    Mu = [];
    Sigma = [];
end

end