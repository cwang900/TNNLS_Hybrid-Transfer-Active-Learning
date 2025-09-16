clear
clc

%%

NDomain = 3; 
NObs = 30; 
NoiseSigma = [0.1 0.03];

Data = load("NC_EIS_Capacity.mat");
Data = Data.Data.T25;
IndState = [1 2 5];

%%

NRep = 100;
NPool = 1000;
NRand = 100;
NTarget = 5;
yTr = cell(1,NDomain);
XTr = cell(1,NDomain);
f = cell(1,NDomain);
ypool = cell(1,NDomain);
Xpool = cell(1,NDomain);
fpool = cell(1,NDomain);
XSource = zeros((NDomain-1)*NObs,NRep);
ySource = zeros((NDomain-1)*NObs,2*NRep);
XTarget = zeros(NTarget,NRep);
yTarget = zeros(NTarget,2*NRep);
XPool = zeros(NPool,NRep);
yPool = zeros(NPool,2*NRep);
yTrue = zeros(NPool,2*NRep);
XRand = zeros(100,NRep);
yRand = zeros(100,2*NRep);

%%

for iRep = 1:NRep
    IndCell = randperm(8,1);
    StateString = {'I','II','III','IV','V','VI','IX'};
    fAll = cell(1,8);
    XAll = cell(1,8);
    ii = 1;
    for iState = IndState
        Temp = ['Cell' num2str(IndCell)];
        State = ['State' StateString{iState}];
        fAll{ii} = [flipud(Data.(Temp).(State).EIS_Re(:,1)) ...
            flipud(Data.(Temp).(State).EIS_Im(:,1))];
        XAll{ii} = flipud(Data.(Temp).(State).Freq(:,1));
        XAll{ii} = repmat(XAll{ii},1,2);
        ii = ii+1;
    end
    X0 = XAll{1};
    Xmin = XAll{1}(1);
    Xmax = XAll{1}(end);
    XXpool = linspace(Xmin,Xmax,NPool)';

    yTr = cell(1,NDomain);
    XTr = cell(1,NDomain);
    f = cell(1,NDomain);
    ypool = cell(1,NDomain);
    Xpool = cell(1,NDomain);
    fpool = cell(1,NDomain);
    for i = 1:NDomain
        if i<NDomain
            p = sort(lhsdesign(NObs,1));
            XX = Xmin+(Xmax-Xmin)*p;
            XTr{i} = XX;
            ff = interp1(XAll{i}(:,1),fAll{i},XX);
            Xpool{i} = [XXpool XXpool];
            fpool{i} = interp1(XAll{i}(:,1),fAll{i},XXpool);
            yTr{i} = ff+normrnd(zeros(NObs,2),NoiseSigma.*ones(NObs,2));
        else
            p = sort(lhsdesign(5,1));
            XX = Xmin+(Xmax-Xmin)*p;
            XTr{i} = XX;
            ff = interp1(XAll{i}(:,1),fAll{i},XX);
            Xpool{i} = XXpool;
            fpool{i} = interp1(XAll{i}(:,1),fAll{i},XXpool);
            yTr{i} = ff+normrnd(zeros(NTarget,2),NoiseSigma.*ones(NTarget,2));
        end
        ypool{i} = fpool{i}+normrnd(zeros(NPool,2),NoiseSigma.*ones(NPool,2));
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
    yTrue(:,(iRep-1)*2+1:iRep*2) = fpool{NDomain};
    
    Xtemp = unifrnd(Xmin,Xmax,NRand,1);
    XRand(:,iRep) = Xtemp;
    yTemp = interp1(XAll{NDomain}(:,1),fAll{NDomain},Xtemp);
    yRand(:,(iRep-1)*2+1:iRep*2) = yTemp+normrnd(zeros(NRand,2),NoiseSigma.*ones(NRand,2));

end

writematrix(XSource,'NN4AL/Data4NN/XSourceCase.xls');
writematrix(ySource,'NN4AL/Data4NN/ySourceCase.xls');
writematrix(XTarget,'NN4AL/Data4NN/XTargetCase.xls');
writematrix(yTarget,'NN4AL/Data4NN/yTargetCase.xls');
writematrix(XPool,'NN4AL/Data4NN/XPoolCase.xls');
writematrix(yPool,'NN4AL/Data4NN/yPoolCase.xls');
writematrix(XRand,'NN4AL/Data4NN/XRandCase.xls');
writematrix(yRand,'NN4AL/Data4NN/yRandCase.xls');
writematrix(yTrue,'NN4AL/Data4NN/yTrueCase.xls');
