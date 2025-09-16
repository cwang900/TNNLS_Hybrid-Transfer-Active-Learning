clear
clc

%%

NDomain = 2; 
NStream = 2;
NObs = 15; 
NTarget = 5;
Xmin = 0; 
Xmax = 5; 

%% 

w0 = 1.5;
    rho = [[0.6; 0.4] [0.4; 0.6] [0.3; 0.7]]; lhsdesign(NStream,NDomain);
    g1 = @(u,x,rho) u+3*((1-rho)*sin(w0*x(:,1))+rho*cos(2*x(:,2)));
    g2 = @(u,x,rho) u+3*((1-rho)*cos(w0*x(:,2))+rho*sin(2*x(:,1)));
    g = @(u,x,eta,rho) u+...
        exp(-(x(:,1)+x(:,2))./eta).*(rho*2*sin(w0*x(:,1))+(1-rho)*2*cos(w0*x(:,2)));

NoiseSigma = 0.1;

NRep = 100;
NPool = 30;

%% 

yTr = cell(1,NDomain);
XTr = cell(1,NDomain);
f = cell(1,NDomain);
ypool = cell(1,NDomain);
Xpool = cell(1,NDomain);
fpool = cell(1,NDomain);
Fun = {g1,g2};
XSource = zeros((NDomain-1)*NObs^2,2*NRep);
ySource = zeros((NDomain-1)*NObs^2,2*NRep);
XTarget = zeros(NTarget^2,2*NRep);
yTarget = zeros(NTarget^2,2*NRep);
XPool = zeros(NPool^2,2*NRep);
yPool = zeros(NPool^2,2*NRep);
XRand = zeros(100,2*NRep);
yRand = zeros(100,2*NRep);
yTrue = zeros(NPool^2,2*NRep);

NObs = 15;
NObsT = 5;

for iRep = 1:NRep
    % 均匀设定x采样点
    tpool = linspace(Xmin,Xmax,NPool)';
    [x1, x2] = meshgrid(tpool,tpool);
    Xpool = [x1(:) x2(:)];
    eta = unifrnd(-10,-5,NStream,1);
    U = unifrnd(5,10,NStream,NDomain);%
    for i = 1:NDomain
        if i<NDomain
            p = sort(lhsdesign(NObs,2));
            [x1, x2] = meshgrid(p(:,1),p(:,2));
            XTr{i} = Xmin+(Xmax-Xmin)*[x1(:) x2(:)];
            f{i} = [Fun{i}(U(1,i),XTr{i},rho(1,i)) Fun{i}(U(2,i),XTr{i},rho(2,i))];
            fpool{i} = [Fun{i}(U(1,i),Xpool,rho(1,i)) Fun{i}(U(2,i),Xpool,rho(2,i))];
        else
            f{i} = [];
            fpool{i} = [];
            p = sort(lhsdesign(NTarget,2));
            [x1, x2] = meshgrid(p(:,1),p(:,2));
            XTr{i} = Xmin+(Xmax-Xmin)*[x1(:) x2(:)];
            for ii = 1:NStream
                f{i} = [f{i} g(U(ii,i),XTr{i},eta(ii),rho(ii,i))];
                fpool{i} = [fpool{i} g(U(ii,i),Xpool,eta(ii),rho(ii,i))];
            end
        end
        yTr{i} = f{i}+normrnd(0,NoiseSigma,size(f{i}));
        ypool = fpool{i}+normrnd(0,NoiseSigma,size(fpool{i}));
    end
    XStemp = XTr{1};%[ XTr{2}];
    ySTemp = yTr{1};%[ yTr{2}];
    XTtemp = XTr{NDomain};
    yTTemp = yTr{NDomain};
    XpTemp = Xpool;
    ypTemp = ypool;
    XSource(:,(iRep-1)*2+1:iRep*2) = XStemp;
    ySource(:,(iRep-1)*2+1:iRep*2) = ySTemp;
    XTarget(:,(iRep-1)*2+1:iRep*2) = XTtemp;
    yTarget(:,(iRep-1)*2+1:iRep*2) = yTTemp;
    XPool(:,(iRep-1)*2+1:iRep*2) = XpTemp;
    yPool(:,(iRep-1)*2+1:iRep*2) = ypTemp;
    yTrue(:,(iRep-1)*2+1:iRep*2) = fpool{NDomain};

    Xtemp = unifrnd(Xmin,Xmax,100,2);
    XRand(:,(iRep-1)*2+1:iRep*2) = Xtemp;
    yTemp = [g(U(1,end),Xtemp,eta(1),rho(1,end)) g(U(2,end),Xtemp,eta(2),rho(2,end))];
    yRand(:,(iRep-1)*2+1:iRep*2) = yTemp+normrnd(0,NoiseSigma,size(yTemp));
end

SourceFolder = [];
writematrix(XSource,[SourceFolder 'XSource2D.xlsx']);
writematrix(ySource,[SourceFolder 'ySource2D.xlsx']);
writematrix(XTarget,[SourceFolder 'XTarget2D.xlsx']);
writematrix(yTarget,[SourceFolder 'yTarget2D.xlsx']);
writematrix(XPool,[SourceFolder 'XPool2D.xlsx']);
writematrix(yPool,[SourceFolder 'yPool2D.xlsx']);
writematrix(XRand,[SourceFolder 'XRand2D.xlsx']);
writematrix(yRand,[SourceFolder 'yRand2D.xlsx']);
writematrix(yTrue,[SourceFolder 'yTrue2D.xlsx']);


