
function K = CovSGPTransfer(theta,X1,X2,SeqX1,SeqX2,NStream)
    D = size(theta,2);
    if D == 1
        theta = repmat(theta,1,2);
    end
    Hyper.Rho = theta(1:2*NStream,:)';
    Hyper.lambda = theta(2*NStream+1:4*NStream,:)';
    K = CovSub(X1,X2,SeqX1,SeqX2,Hyper,NStream);

    [m,n] = size(K);
    if D == 1 && m == n && sum((X1-X2).^2,"all") == 0
        sigma2 = theta(end-NStream+1:end,1).^2';
        sigma2 = repelem(sigma2,SeqX1);
        sigma2 = diag(sigma2);
        K = K+sigma2;
    end

end

function K = CovSub(X1,X2,SeqX1,SeqX2,Hyper,NStream)

rho_ij = Hyper.Rho; % [rho(d1)_i0 rho(d1)_ii rho(d1)_i'0 rho(d1)_i'i';
                    %  rho(d2)_i0 rho(d2)_ii rho(d2)_i'0 rho(d2)_i'i']
                    %  d1/d2--domain index
lambda_ij = (Hyper.lambda).^2; % (exp(Hyper.lambda)).^2; % [l(d1)_i0 l(d1)_ii l(d1)_i'0 l(d1)_i'i';
                          %  l(d2)_i0 l(d2)_ii l(d2)_i'0 l(d2)_i'i']
                          %  d1/d2--domain index
beta = 0.5^2;

%% Prepare scale parameter

% lambda_square_i0j0 = sqrt(lambda_ij(1,[1 3]).^2'+lambda_ij(2,[1 3]).^2);
% % lambda_square_i0j0 = diag(diag(lambda_square_i0j0));
% 
% scale_const_i0j0 = sqrt(2*(lambda_ij(1,[1 3])'.*lambda_ij(2,[1 3])))./...
%                                                 lambda_square_i0j0;

% lambda_square_i0j0 = sqrt(lambda_ij(1,[1 3])'+lambda_ij(2,[1 3])+beta);
% 
% scale_const_i0j0 = sqrt(2)*((lambda_ij(1,[1 3])'.*lambda_ij(2,[1 3]).*beta).^(1/4))./...
%                                                 lambda_square_i0j0;

lambda_square_i0j0 = sqrt(lambda_ij(1,1:2:NStream*2)'+lambda_ij(2,1:2:NStream*2));
scale_const_i0j0 = sqrt(2)*((lambda_ij(1,1:2:NStream*2)'.*lambda_ij(2,1:2:NStream*2))).^(1/4)./...
                                                (lambda_square_i0j0);

scale_const_i0j0 = diag(diag(scale_const_i0j0));

lambda_square_i0j0 = repelem(lambda_square_i0j0,SeqX1,SeqX2);
scale_const_i0j0 = repelem(scale_const_i0j0,SeqX1,SeqX2);

if rho_ij(1,2) ~= 0
    lambda_square_iijj = sqrt(lambda_ij(1,2:2:NStream*2).^2'+lambda_ij(2,2:2:NStream*2).^2);

    scale_const_iijj = sqrt(2*(lambda_ij(1,2:2:NStream*2)'.*lambda_ij(2,2:2:NStream*2)))./...
        lambda_square_iijj;

    lambda_square_iijj = repelem(lambda_square_iijj,SeqX1,SeqX2);
    scale_const_iijj = repelem(scale_const_iijj,SeqX1,SeqX2);
end


%% Prepare correlation parameter

rho_square_i0j0 = rho_ij(1,1:2:NStream*2)'.*rho_ij(2,1:2:NStream*2);
rho_square_i0j0 = diag(diag(rho_square_i0j0));
rho_square_i0j0 = repelem(rho_square_i0j0,SeqX1,SeqX2);

if rho_ij(1,2) ~= 0
    rho_square_iijj = rho_ij(1,2:2:NStream*2)'.*rho_ij(2,2:2:NStream*2);
    rho_square_iijj(1,2) = 0;
    rho_square_iijj(2,1) = 0;
    rho_square_iijj = repelem(rho_square_iijj,SeqX1,SeqX2);
end

%% Prepare the input

D2XX = sq_dist(X1', X2');

Ki0j0 = rho_square_i0j0.*scale_const_i0j0.*exp(-1/2*D2XX./lambda_square_i0j0.^2);
if rho_ij(1,2) ~= 0
    Kiijj = rho_square_iijj.*scale_const_iijj.*exp(-1/2*D2XX./lambda_square_iijj.^2);
    K = Ki0j0+Kiijj;
else
    K = Ki0j0;
end


end






