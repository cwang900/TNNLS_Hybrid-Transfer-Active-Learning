function theta = OptJOffline2D(theta0,y,Xf,XInd,CovType,NStream)

opt.algorithm = NLOPT_LD_MMA;%NLOPT_LN_AUGLAG_EQ;
ub = inf*ones(1,length(theta0));
lb = -inf*ones(1,length(theta0));
NDomain = length(y);
% if strcmp(CovType,'CovOneShare')
%     lb(4*(NDomain-1)+2*NDomain+(1:2*NDomain)) = 0.1;
% end
lb(end) = 0.05;
ub(end) = 0.5;
opt.lower_bounds = lb;
opt.upper_bounds = ub;
opt.xtol_rel = 1e-4;

% options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);

IndFix = [];
ParFix = [];
NIter = 100+round((length(y{end})-50)^1.5*10);
opt.maxeval = 5000;%NIter;%
% if strcmp(CovType,'CovOneShare')
%     opt.maxeval = 3000;
% else
%     opt.maxeval = 1000;
% end
opt.min_objective = @(theta) ObjSubFix(theta,y,Xf,XInd,CovType,ParFix,IndFix,NStream);
theta = nlopt_optimize(opt, theta0);
    
end

%%

function [J,dJ] = ObjSubFix(theta,y,Xf,XInd,CovType,ParFix,IndFix,NStream)

% theta0 = theta;
% tic
[J,dJ] = ObjSub(theta,y,Xf,XInd,CovType,ParFix,IndFix,NStream);
% toc

end

%%

function [J,dJ] = ObjSub(theta,y,Xf,XInd,CovType,ParFix,IndFix,NStream)

theta0 = theta;
% theta = zeros(1,length(IndFix));
% theta(IndFix) = theta0;
% theta(~IndFix) = ParFix;
NDomain = length(y);

%%

Aim = 'Optimization';
if NDomain == 1
    [~, Inv_Cxx_nn, ~, Det] = ...
        CovCalAssembled2D(theta,y,Xf,CovType,Aim,NStream);
    yKL = Inv_Cxx_nn*y{1};
    J = 1/2*Det+1/2*(yKL'*yKL);
else
    [Myn, Inv_Cxx_nn, Inv_KS, Det] = ...
        CovCalAssembled2D(theta,y,Xf,CovType,Aim,NStream);

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

