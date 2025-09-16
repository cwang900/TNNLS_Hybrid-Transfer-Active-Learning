clear
clc

%%

NDomain = 8; 
NObs = 30; %
NTarget = 5;
Xmin = 0; 
Xmax = 5; 

X0 = linspace(Xmin,Xmax,1000)';
y = cell(1,NDomain);
X = cell(1,NDomain);

%% Parametric simulation


g1 = @(b,x) b(1)+b(2)*x+b(3)*x.^2+b(4)*x.^3;
g2 = @(w,rho,x) 2*(rho*sin(w*x)+(1-rho)*cos(w*x+pi/3));
g = @(rho,x) 0*x.^3+1*x.^2-x+2+4*(rho*sin(1*x)+(1-rho)*cos(1*x));

Fun = {g1,g2};

NoiseSigma = 1;

NRep = 100;
NPool = 1000;

XSource = zeros((NDomain-1)*NObs,NRep);
ySource = zeros((NDomain-1)*NObs,2*NRep);
XTarget = zeros(NTarget,NRep);
yTarget = zeros(NTarget,2*NRep);
XPool = zeros(NPool,NRep);
yPool = zeros(NPool,2*NRep);
XRand = zeros(100,NRep);
yRand = zeros(100,2*NRep);
yTrue = zeros(NPool,2*NRep);

for iRep = 1:NRep

    b = [unifrnd(0,5,1,NDomain); unifrnd(0,2,3,NDomain);];%
    b(4,:) = 0;
    Ind = randi([1 4],1,NDomain);
    Ind = sub2ind(size(b),Ind,1:NDomain);
    b(Ind) = 0;
    w = normrnd(1,0.1,1,NDomain);
    rho = [1 0; 0.2 0.8];%

    %
    p = sort(lhsdesign(NObs,NDomain));
    t = Xmin+(Xmax-Xmin)*p;
    t = mat2cell(t,NObs,ones(1,NDomain));
    p = sort(lhsdesign(5,1));
    t{end} = Xmin+(Xmax-Xmin)*p;
    tpool = linspace(Xmin,Xmax,NPool)';

    yTr = cell(1,NDomain);
    XTr = cell(1,NDomain);
    f = cell(1,NDomain);
    ypool = cell(1,NDomain);
    Xpool = cell(1,NDomain);
    fpool = cell(1,NDomain);
    for i = 1:NDomain
        XTr{i} = t{i};
        Xpool{i} = tpool;
        if i<NDomain
            f{i} = [g1(b(:,i),t{i})+g2(w(i),rho(1,1),t{i}) ...
                g1([-b(1,i); b(2:end,i)],t{i})+g2(w(i),rho(1,2),t{i})];
            fpool{i} = [g1(b(:,i),tpool)+g2(w(i),rho(1,1),tpool) ...
                g1([-b(1,i); b(2:end,i)],tpool)+g2(w(i),rho(1,2),tpool)];
        else
            f{i} = [g(rho(2,1),t{i})  g(rho(2,2),t{i})];
            fpool{i} = [g(rho(2,1),tpool)  g(rho(2,2),tpool)];
        end
        yTr{i} = f{i}+normrnd(0,NoiseSigma,size(f{i}));
        ypool{i} = fpool{i}+normrnd(0,NoiseSigma,size(fpool{i}));
    end
    XStemp = [XTr{1}; XTr{2}; XTr{3}; XTr{4}; XTr{5}; XTr{6}; XTr{7}];
    ySTemp = [yTr{1}; yTr{2}; yTr{3}; yTr{4}; yTr{5}; yTr{6}; yTr{7}];
    XTtemp = XTr{NDomain};
    yTTemp = yTr{NDomain};
    XpTemp = Xpool{NDomain};
    ypTemp = ypool{NDomain};
    XSource(:,iRep) = XStemp;
    ySource(:,(iRep-1)*2+1:iRep*2) = ySTemp;
    XTarget(:,iRep) = XTtemp;
    yTarget(:,(iRep-1)*2+1:iRep*2) = yTTemp;
    XPool(:,iRep) = XpTemp;
    yPool(:,(iRep-1)*2+1:iRep*2) = ypTemp;
    yTrue(:,(iRep-1)*2+1:iRep*2) = fpool{NDomain};

    Xtemp = unifrnd(Xmin,Xmax,100,1);
    XRand(:,iRep) = Xtemp;
    yTemp = [g(rho(2,1),Xtemp) g(rho(2,2),Xtemp)];
    yRand(:,(iRep-1)*2+1:iRep*2) = yTemp+normrnd(0,NoiseSigma,size(yTemp));

end
writematrix(XSource,'NN4AL/Data4NN/XSource8Domain.xls');
writematrix(ySource,'NN4AL/Data4NN/ySource8Domain.xls');
writematrix(XTarget,'NN4AL/Data4NN/XTarget8Domain.xls');
writematrix(yTarget,'NN4AL/Data4NN/yTarget8Domain.xls');
writematrix(XPool,'NN4AL/Data4NN/XPool8Domain.xls');
writematrix(yPool,'NN4AL/Data4NN/yPool8Domain.xls');
writematrix(XRand,'NN4AL/Data4NN/XRand8Domain.xls');
writematrix(yRand,'NN4AL/Data4NN/yRand8Domain.xls');
writematrix(yTrue,'NN4AL/Data4NN/yTrue8Domain.xls');
