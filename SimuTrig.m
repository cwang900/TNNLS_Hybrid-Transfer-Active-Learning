clear
clc

%%

NDomain = 3; 
NObs = 30; 
NTarget = 5;
Xmin = -5; 
Xmax = 5; 

%% 

w0 = 0.7;
g1 = @(u,x,rho) u+2*((1-rho)*sin(w0*x)+rho*sin(1*x));
g2 = @(u,x,rho) u+2*((1-rho)*cos(w0*x)+rho*cos(1*x));
g = @(u,x,Phi,eta,rho) u+...
    exp(-x./eta).*(rho*2*sin(w0*x)+(1-rho)*2*cos(w0*x));

rho = [[0.5; 0.7] [0.5; 0.7] [0.7; 0.3]];%[0.8 0.2 0.6];%
eta = [5,-10];
NoiseSigma = 0.5;

NRep = 100;
NPool = 1000;

%% 生成数据

yTr = cell(1,NDomain);
XTr = cell(1,NDomain);
f = cell(1,NDomain);
ypool = cell(1,NDomain);
Xpool = cell(1,NDomain);
fpool = cell(1,NDomain);
Fun = {g1,g2};
XSource = zeros((NDomain-1)*NObs,NRep);
ySource = zeros((NDomain-1)*NObs,2*NRep);
XTarget = zeros(NTarget,NRep);
yTarget = zeros(NTarget,2*NRep);
XPool = zeros(NPool,NRep);
yPool = zeros(NPool,2*NRep);
XRand = zeros(100,NRep);
yRand = zeros(100,2*NRep);
yTrue = zeros(NPool,2*NRep);
XX = zeros(50,NRep);
yy = zeros(50,2*NRep);

for iRep = 1:NRep
    % 均匀设定x采样点
    p = sort(lhsdesign(NObs,NDomain));
    t = Xmin+(Xmax-Xmin)*p;
    t = mat2cell(t,NObs,ones(1,NDomain));
    p = sort(lhsdesign(NTarget,1));
    t{end} = Xmin+(Xmax-Xmin)*p;
    tpool = linspace(Xmin,Xmax,NPool)';
    U = [unifrnd(0,5,1,NDomain);unifrnd(-5,0,1,NDomain);];
    Phi = [0,pi/4];%normrnd(2,1,1,2);
    for i = 1:NDomain
        XTr{i} = t{i};
        Xpool{i} = tpool;
        if i<NDomain
            f{i} = [Fun{i}(U(1,i),t{i},rho(1,i)) Fun{i}(U(2,i),t{i},rho(2,i))];
            fpool{i} = [Fun{i}(U(1,i),tpool,rho(1,i)) Fun{i}(U(2,i),tpool,rho(2,i))];
        else
            f{i} = [g(U(1,i),t{i},Phi(1),eta(1),rho(1,i)) g(U(2,i),t{i},Phi(2),eta(2),rho(2,i))];
            fpool{i} = [g(U(1,i),tpool,Phi(1),eta(1),rho(1,i)) g(U(2,i),tpool,Phi(2),eta(2),rho(2,i))];
        end
        yTr{i} = f{i}+normrnd(0,NoiseSigma,size(f{i}));
        ypool{i} = fpool{i}+normrnd(0,NoiseSigma,size(fpool{i}));
    end
    XStemp = [XTr{1}; XTr{2}];
    ySTemp = [yTr{1}; yTr{2}];
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
    Xtemp = unifrnd(Xmin,Xmax,100,1);
    yTrue(:,(iRep-1)*2+1:iRep*2) = fpool{NDomain};

    XRand(:,iRep) = Xtemp;
    yTemp = [g(U(1,3),Xtemp,Phi(1),eta(1),rho(1,3)) g(U(2,3),Xtemp,Phi(2),eta(2),rho(2,3))];
    yRand(:,(iRep-1)*2+1:iRep*2) = yTemp+normrnd(0,NoiseSigma,size(yTemp));
end

writematrix(XSource,'NN4AL/Data4NN/XSource.xls');
writematrix(ySource,'NN4AL/Data4NN/ySource.xls');
writematrix(XTarget,'NN4AL/Data4NN/XTarget.xls');
writematrix(yTarget,'NN4AL/Data4NN/yTarget.xls');
writematrix(XPool,'NN4AL/Data4NN/XPool.xls');
writematrix(yPool,'NN4AL/Data4NN/yPool.xls');
writematrix(XRand,'NN4AL/Data4NN/XRand.xls');
writematrix(yRand,'NN4AL/Data4NN/yRand.xls');
writematrix(yTrue,'NN4AL/Data4NN/yTrue.xls');


