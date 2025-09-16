function theta = OptJOffline(theta0,y,Xf,XInd,CovType,NStream)

opt.algorithm = NLOPT_LD_MMA;%NLOPT_LN_AUGLAG_EQ;
ub = inf*ones(1,length(theta0));
lb = -inf*ones(1,length(theta0));
NDomain = length(y);

lb(end) = 0.1;
ub(end) = 1;
opt.lower_bounds = lb;
opt.upper_bounds = ub;
opt.xtol_rel = 1e-3;

% options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);

IndFix = [];
ParFix = [];
opt.maxeval = 10000;

opt.min_objective = @(theta) ObjSubFix(theta,y,Xf,XInd,CovType,ParFix,IndFix,NStream);
theta = nlopt_optimize(opt, theta0);
    
end

%%

function [J,dJ] = ObjSource(theta,y,Xf,Xu,CovType,ParFix)

theta0 = theta;
theta = [theta ParFix];
NDomain = length(y);

%% Prepare model parameters

rhoSource = reshape(theta(1:2*(NDomain)),2,[]);
rhoSource = repelem(rhoSource,2,1);
rhoSource([2 4],:) = 0;
theta(1:2*(NDomain)) = [];

lambdaSource = reshape(theta(1:2*(NDomain)),2,[]);
lambdaSource = repelem(lambdaSource,2,1);
theta(1:2*(NDomain)) = [];

sigmaSource = repmat(theta(1:1*(NDomain)),2,1);

thetaSource = [rhoSource;lambdaSource;sigmaSource];
thetaSource = mat2cell(thetaSource,10,ones(1,NDomain));

%% Covariance calculation

CovFun = str2func(CovType);

XSource = Xf;

SeqXSource = cellfun(@length,XSource)/2;
SeqXSource = repmat(SeqXSource,2,1);
SeqXSource = mat2cell(SeqXSource,2,ones(1,NDomain));

KSource = cellfun(CovFun,thetaSource,XSource,XSource,...
        SeqXSource,SeqXSource,'UniformOutput',false);
KL = cellfun(@chol,KSource,'UniformOutput',false);
Det = @(K) 2*sum(log(diag(K)));
DetSource = cellfun(Det,KL);
DetSource = sum(DetSource);

Inv_KSL = cellfun(@inv,KL,'UniformOutput',false);
Fun = @(K,y) (K'*y)'*(K'*y);
yKSy = cellfun(Fun,Inv_KSL,y);

J = 1/2*DetSource+1/2*sum(yKSy);

%% Calculate objective function

if nargout == 2

    fun = @(theta)ObjSource(theta,y,Xf,Xu,CovType,ParFix);
    dJ = Grad(fun, theta0);

end

end

%%

function [J,dJ] = ObjSubFix(theta,y,Xf,XInd,CovType,ParFix,IndFix,NStream)

% theta0 = theta;
[J,dJ] = ObjSub(theta,y,Xf,XInd,CovType,ParFix,IndFix,NStream);

end

%%

function [J,dJ] = ObjSub(theta,y,Xf,XInd,CovType,ParFix,IndFix,NStream)

theta0 = theta;
% theta = zeros(1,length(IndFix));
% theta(IndFix) = theta0;
% theta(~IndFix) = ParFix;
NDomain = length(y);
if strcmp(CovType,'CovLMC')
    if length(theta) == (2*NStream+2)*NDomain-NStream
        rhoSource = theta(1:NStream*(NDomain-1));
        theta(1:NStream*(NDomain-1)) = [];

        rhoTarget = theta(1:NStream*(NDomain));
        theta(1:NStream*(NDomain)) = [];

        lambdaSource = theta(1:NDomain-1);
        lambdaSource = repelem(lambdaSource,NStream,1);
        theta(1:NDomain-1) = [];

        lambdaTarget = theta(1);
        lambdaTarget = repelem(lambdaTarget,NStream,1);
        lambdaTarget = [lambdaSource lambdaTarget];
        theta(1) = [];

        sigma = theta;
    else
        rhoSource = theta(1:NStream*(NDomain-1));
        theta(1:NStream*(NDomain-1)) = [];

        rhoTarget = theta(1:NStream*(NDomain+1));
        theta(1:NStream*(NDomain+1)) = [];

        lambdaSource = theta(1:NDomain-1);
        lambdaSource = repelem(lambdaSource,NStream,1);
        theta(1:NDomain-1) = [];

        lambdaTarget = theta(1:2)';
        lambdaTarget = repmat(lambdaTarget,NStream/2,1);
        lambdaTarget = [lambdaSource lambdaTarget];
        theta(1:2) = [];

        sigma = theta;
    end
    
    theta = [rhoSource rhoTarget lambdaSource(:)' lambdaTarget(:)' sigma];
end

%%

Aim = 'Optimization';
if NDomain == 1
    [~, Inv_Cxx_nn, ~, Det] = ...
        CovCalAssembledVV(theta,y,Xf,CovType,Aim,NStream);
    yKL = Inv_Cxx_nn*y{1};
    J = 1/2*Det+1/2*(yKL'*yKL);
else
    [Myn, Inv_Cxx_nn, Inv_KS, Det] = ...
        CovCalAssembledVV(theta,y,Xf,CovType,Aim,NStream);

    yS = y(1:NDomain-1);
    Prod = @(K,y) y'*K*y;
    yKSy = cellfun(Prod,Inv_KS,yS);
    yKSy = sum(yKSy);
    
    yKT = Inv_Cxx_nn*(y{end}-Myn);
    yKTy = yKT'*yKT;

    J = 1/2*Det+1/2*(yKSy+yKTy);
end


%% Calculate objective function

if nargout == 2

    fun = @(theta)ObjSub(theta,y,Xf,XInd,CovType,ParFix,IndFix,NStream);
    dJ = Grad(fun, theta0);

end

end

