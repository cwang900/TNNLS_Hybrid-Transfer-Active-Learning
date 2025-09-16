
function K = CovLMC(theta,X1,X2,SeqX1,SeqX2,NStream)
    D = size(theta,2);
    if D == 1
        theta = repmat(theta,1,2);
    end
    Hyper.Rho = theta(1:2*NStream,:)';
    Hyper.lambda = theta(2*NStream+1:4*NStream,:);
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

D = size(X1,2); % input dimension
% rho_mm = exp(Hyper.Rho);
rho_mm = Hyper.Rho;

if rho_mm(1,2) == 0
    scale_mm = abs(Hyper.lambda);
    scale_mm = scale_mm(1:2:2*NStream,1);
    scale_mm = reshape(scale_mm,D,[]);

    X1scale_mm = repelem(scale_mm',SeqX1,1);
    X2scale_mm = repelem(scale_mm',SeqX2,1);
    X1_mm = X1./X1scale_mm/sqrt(2);
    X2_mm = X2./X2scale_mm/sqrt(2);

    D2mm = sq_dist(X1_mm', X2_mm');
    rho_mm_square = rho_mm(1,1:2:2*NStream)'.*rho_mm(2,1:2:2*NStream);
    % rho_mm_square = diag(diag(rho_mm_square));
    rho_mm_square = repelem(rho_mm_square,SeqX1,SeqX2);

    K = rho_mm_square.*exp(-1/2*D2mm);
else
    K = 0;
    for i = 1:2
        scale_mm = abs(Hyper.lambda);
        scale_mm = scale_mm(i:2:2*NStream,1);
        scale_mm = reshape(scale_mm,D,[]);

        X1scale_mm = repelem(scale_mm',SeqX1,1);
        X2scale_mm = repelem(scale_mm',SeqX2,1);
        X1_mm = X1./X1scale_mm/sqrt(2);
        X2_mm = X2./X2scale_mm/sqrt(2);

        D2mm = sq_dist(X1_mm', X2_mm');
        rho_mm_square = rho_mm(1,i:2:2*NStream)'.*rho_mm(2,i:2:2*NStream);
        % rho_mm_square = diag(diag(rho_mm_square));
        rho_mm_square = repelem(rho_mm_square,SeqX1,SeqX2);

        K = K+rho_mm_square.*exp(-1/2*D2mm);
    end
end


end






