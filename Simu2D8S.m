clear
clc

%%

NDomain = 2; 
NObs = 30; 
NTarget = 5;
Xmin = -5; 
Xmax = 5; 
NStream = 8;

%% 

w0 = 0.7;
g1 = @(u,x,rho) u+5*((1-rho)*sin(w0*x)+rho*sin(1*x));
g2 = @(u,x,rho) u+5*((1-rho)*cos(w0*x)+rho*sin(1*x));
g = @(u,x,eta,rho) u+...
    exp(-x./eta).*(rho*2*sin(w0*x)+(1-rho)*2*cos(w0*x));

NoiseSigma = 0.5;

NRep = 100;
NPool = 1000;

%% 

yTr = cell(1,NDomain);
XTr = cell(1,NDomain);
f = cell(1,NDomain);
ypool = cell(1,NDomain);
Xpool = cell(1,NDomain);
fpool = cell(1,NDomain);
Fun = {g1,g2};
XSource = zeros((NDomain-1)*NObs,NRep);
ySource = zeros((NDomain-1)*NObs,NStream*NRep);
XTarget = zeros(NTarget,NRep);
yTarget = zeros(NTarget,NStream*NRep);
XPool = zeros(NPool,NRep);
yPool = zeros(NPool,NStream*NRep);
XRand = zeros(100,NRep);
yRand = zeros(100,NStream*NRep);
yTrue = zeros(NPool,NStream*NRep);
XX = zeros(50,NRep);
yy = zeros(50,NStream*NRep);

for iRep = 1:NRep
    
    p = sort(lhsdesign(NObs,NDomain));
    t = Xmin+(Xmax-Xmin)*p;
    t = mat2cell(t,NObs,ones(1,NDomain));
    p = sort(lhsdesign(NTarget,1));
    t{end} = Xmin+(Xmax-Xmin)*p;
    tpool = linspace(Xmin,Xmax,NPool)';
    rho = lhsdesign(NStream,NDomain);
    eta = [unifrnd(-10,-5,NStream,1);unifrnd(4,5,NStream,1)];%[5,-10];
    U = unifrnd(-5,5,NStream,NDomain);%zeros(2,NDomain);%
    for i = 1:NDomain
        f{i} = [];
        fpool{i} = [];
        XTr{i} = t{i};
        Xpool{i} = tpool;
        if i<NDomain
            for ii = 1:NStream
                f{i} = [f{i} g1(U(ii,i),t{i},rho(ii,i))];
                fpool{i} = [fpool{i} g1(U(ii,i),tpool,rho(ii,i))];
            end
        else
            for ii = 1:NStream
                f{i} = [f{i} g(U(ii,i),t{i},eta(ii),rho(ii,i))];
                fpool{i} = [fpool{i} g(U(ii,i),tpool,eta(ii),rho(ii,i))];
            end
        end
        yTr{i} = f{i}+normrnd(0,NoiseSigma,size(f{i}));
        ypool{i} = fpool{i}+normrnd(0,NoiseSigma,size(fpool{i}));
    end
    XStemp = cell2mat(XTr(1:NDomain-1));
    ySTemp = cell2mat(yTr(1:NDomain-1));
    XTtemp = XTr{NDomain};
    yTTemp = yTr{NDomain};
    XpTemp = Xpool{NDomain};
    ypTemp = ypool{NDomain};
    XSource(:,iRep) = XStemp;
    ySource(:,(iRep-1)*NStream+1:iRep*NStream) = ySTemp;
    XTarget(:,iRep) = XTtemp;
    yTarget(:,(iRep-1)*NStream+1:iRep*NStream) = yTTemp;
    XPool(:,iRep) = XpTemp;
    yPool(:,(iRep-1)*NStream+1:iRep*NStream) = ypTemp;
    Xtemp = unifrnd(Xmin,Xmax,100,1);
    yTrue(:,(iRep-1)*NStream+1:iRep*NStream) = fpool{NDomain};

    XRand(:,iRep) = Xtemp;
    yTemp = [];
    for ii = 1:NStream
        yTemp = [yTemp g(U(ii,NDomain),Xtemp,eta(ii),rho(ii,NDomain))];
    end
    yRand(:,(iRep-1)*NStream+1:iRep*NStream) = yTemp+normrnd(0,NoiseSigma,size(yTemp));
end

% set(gcbf,'MaxColumns',1000)
writematrix(XSource,'NN4AL/Data4NN/XSource8Streams.xlsx');
writematrix(ySource,'NN4AL/Data4NN/ySource8Streams.xlsx');
writematrix(XTarget,'NN4AL/Data4NN/XTarget8Streams.xlsx');
writematrix(yTarget,'NN4AL/Data4NN/yTarget8Streams.xlsx');
writematrix(XPool,'NN4AL/Data4NN/XPool8Streams.xlsx');
writematrix(yPool,'NN4AL/Data4NN/yPool8Streams.xlsx');
writematrix(XRand,'NN4AL/Data4NN/XRand8Streams.xlsx');
writematrix(yRand,'NN4AL/Data4NN/yRand8Streams.xlsx');
writematrix(yTrue,'NN4AL/Data4NN/yTrue8Streams.xlsx');


