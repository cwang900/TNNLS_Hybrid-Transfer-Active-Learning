%% Marginal particle filter
function [MuTheta,VarTheta,Index,Muh_t,Ph_t,Estimate] = ...
    MargPartFilter(yt,Xt,XtInd,Particle,Muh_t_1,Ph_t_1,NDomain,CovType,NStream)

NParticle = size(Particle,1);
WNew = zeros(NParticle,1);
Muht = zeros(length(yt),NParticle);
Pht = zeros(length(yt),length(yt),NParticle);
Muh_t = zeros(size(Muh_t_1));
Ph_t = zeros(size(Ph_t_1));

%% Importance sampling

for i = 1:NParticle
    [Muht_Pred, Pht_Pred, Muh_t(:,i), Ph_t(:,:,i)] = PostInd(Particle(i,:)',{yt},{Xt},{Xt}, ...
            Muh_t_1(:,i),Ph_t_1(:,:,i),{XtInd},NDomain,CovType,NStream);
    Lamda = Particle(i,end)^2*eye(length(Muht_Pred));
    WNew(i) = LogLik(yt,Muht_Pred,Pht_Pred,Lamda); % Calculate the weight for each particle
end

% Normalize the weight
% WNew(isnan(WNew)) = -1000;
WNewTemp = exp(WNew);
for i = 1:length(WNewTemp)
    if WNewTemp(i)>10^8 || isnan(WNewTemp(i))
        WNewTemp(i) = unifrnd(10^5,10^7,1,1);
    end
    if WNewTemp(i)<1e-5
        WNewTemp(i) = unifrnd(1e-5,1e-4,1,1);
    end
end
% WNew = WNewTemp;
if ~isnan(sum(WNewTemp))
    WNew = WNewTemp;
% else
%     for i = 1:length(WNewTemp)
%         if WNewTemp(i)>10^4 || isnan(WNewTemp(i))
%             WNewTemp(i) = unifrnd(10^3,10^4,1,1);
%         end
%     end
%     WNew = WNewTemp;
end
% if sum(WNew<0)
%     aa;
% end
WNew = WNew./sum(WNew);

% if isnan(sum(WNew))
%     aa;
% end

%% Weighted estimation

Estimate.Muh_t_hat = Muh_t*WNew;
% Estimate.Muht_hat = Muht*WNew;
MuTheta = Particle'*WNew;
VarTheta = (Particle-MuTheta')'*(WNew.*(Particle-MuTheta'));
Estimate.Theta_hat = MuTheta;
WNew_temp = reshape(WNew,1,1,[]);
Estimate.Ph_t_hat = sum(Ph_t.*WNew_temp,3);
% Estimate.Pht_hat = sum(Pht.*WNew_temp,3);

%% Resampling

% Index = randsample(NParticle,NParticle,true,WNew);
try 
    Index = randsample(NParticle,NParticle,true,WNew(:));
catch ME
    disp(ME.message);
end
end

function lik = LogLik(yt,Mut,Ktt,Lamda)
    Sigma = Ktt+Lamda;
    CholS = chol(Sigma,'lower');
    yLS = CholS\(yt-Mut);
    
    LogDet = 2*sum(log(diag(CholS)));
    lik = -length(yt)/2*log(2*pi)-1/2*LogDet-...
        1/2*(yLS'*yLS);
end
% Particle = Particle(Index,:);

