import os
import torch
import pandas as pd
import torch.nn as nn
import numpy as np
import time
import matplotlib.pyplot as plt
from bmdal_reg.bmdal.feature_data import TensorFeatureData
from bmdal_reg.bmdal.algorithms import select_batch

os.chdir(r'...')   # 将笔记本中的工作目录更改为主目录

xS = pd.read_excel('.../XSource2D.xlsx', header=None)
yS = pd.read_excel('.../ySource2D.xlsx', header=None)
xT = pd.read_excel('.../XTarget2D.xlsx', header=None)
yT = pd.read_excel('.../yTarget2D.xlsx', header=None)
xP = pd.read_excel('.../XPool2D.xlsx', header=None)
yP = pd.read_excel('.../yPool2D.xlsx', header=None)
xR = pd.read_excel('.../XRand2D.xlsx', header=None)
yR = pd.read_excel('.../yRand2D.xlsx', header=None)
yTrue = pd.read_excel('.../yTrue2D.xlsx', header=None)

m=20
class CustomModel(nn.Module):
    def __init__(self):
        super(CustomModel, self).__init__()
        self.layer1 = nn.Linear(2, m)  # 第一层
        self.act1 = nn.SiLU()  # 第一层激活
        self.layer2 = nn.Linear(m, m)  # 第二层
        self.act2 = nn.SiLU()  # 第二层激活
        self.layer3 = nn.Linear(m, 1)  # 输出层1
        self.layer4 = nn.Linear(m, 1)  # 输出层2

    def forward(self, x):
        x = self.layer1(x)  # 第一层线性变换
        x = self.act1(x)     # 激活函数
        x = self.layer2(x)   # 第二层线性变换
        x = self.act2(x)     # 激活函数
        x1 = self.layer3(x)   # 输出1
        x2 = self.layer4(x)  # 输出2
        x = torch.cat((x1, x2), dim=1)  # 使用torch.cat直接将x1和x2合并成一个2列的张量
        return x

NRep = 100 # Number of replications
NRound = 30 # Number of learning round
NSPretrain = 5000 # Number of training steps with source data
NTOffline = 1000 # Number of training steps with offline target data
NTOnline = 500 # Number of training steps during online learning
NStream = 2 # Number of Streams
wl_loss = np.zeros((NRound, NRep))
wl_mape = np.zeros((NRound, NRep))
ExeTime = np.zeros((NRound, NRep))

# SelectionMethod = 1 # for random selection
SelectionMethod = 2  # for active selection
for l in range(NRep):
    # NTOnline = 100
    # Pretrain with source data

    ModelPretrain = CustomModel() # load model
    x_train = torch.tensor(xS.values).float()[:,NStream*l:NStream*l+NStream]#.unsqueeze(1) # load data
    y_train = torch.tensor(yS.values).float()[:,NStream*l:NStream*l+NStream]

    opt = torch.optim.AdamW(ModelPretrain.parameters(), lr=2e-2)
    for epoch in range(NSPretrain):
        y_pred = ModelPretrain(x_train)
        loss1 = ((y_pred[:, 0] - y_train[0:, 0]) ** 2).mean()
        loss2 = ((y_pred[:, 1] - y_train[0:, 1]) ** 2).mean()
        loss = (loss1 + loss2) / 2
        train_rmse = loss.sqrt().item()
        loss.backward()
        opt.step()
        opt.zero_grad()

    # model_filename = f'D:/OneDrive - USTC/Research/Chao/05_Active learning/Code/NN4AL/ModelSource{l}.pth'
    # torch.save(ModelPretrain.state_dict(), model_filename) # save the pretrained model

    # Training with offline target data

    x_train = torch.tensor(xT.values).float()[:,NStream*l:NStream*l+NStream]#.unsqueeze(1)
    y_train = torch.tensor(yT.values).float()[:,NStream*l:NStream*l+NStream]
    x_pool = torch.tensor(xP.values).float()[:,NStream*l:NStream*l+NStream]#.unsqueeze(1)
    y_pool = torch.tensor(yP.values).float()[:,NStream*l:NStream*l+NStream]
    y_true = torch.tensor(yTrue.values).float()[:,NStream*l:NStream*l+NStream]
    xx = torch.tensor(xP.values).float()[::20,NStream*l:NStream*l+NStream]#.unsqueeze(1)
    # yy = torch.tensor(yP.values).float()[::100,NStream*l:NStream*l+NStream]
    # print(xx)
    # print(x_train)

    # print(yP)
    # yp_pred = ModelPretrain(x_pool)
    # plt.plot(x_pool[:, 0].numpy(), y_pool[:, 0].numpy(), '.', color='#BBBBBB')
    # plt.plot(x_pool[:, 0].numpy(), y_pool[:, 1].numpy(), '.', color='#BBBBBB')
    # plt.plot(x_pool[:, 0].numpy(), yp_pred[:, 0].detach().numpy(), '.', color='k')
    # plt.plot(x_pool[:, 0].numpy(), yp_pred[:, 1].detach().numpy(), '.', color='k')
    # plt.show()

    # ModelTarget = CustomModel()
    # model_filename = f'D:/OneDrive - USTC/Research/Chao/05_Active learning/Code/NN4AL/ModelSource{l}.pth'
    # ModelTarget.load_state_dict(torch.load(model_filename))
    tTemp = []
    start = time.time()
    ModelTarget = ModelPretrain
    opt = torch.optim.AdamW(ModelTarget.parameters(), lr=2e-2)
    for epoch in range(NTOffline):
        y_pred = ModelTarget(x_train)
        loss1 = ((y_pred[:, 0] - y_train[0:, 0]) ** 2).mean()
        loss2 = ((y_pred[:, 1] - y_train[0:, 1]) ** 2).mean()
        loss = (loss1 + loss2) / 2
        loss.backward()
        opt.step()
        opt.zero_grad()

    train_rmse = loss.sqrt().item()
    yp_pred = ModelTarget(x_pool)
    # print(y_pool)
    pool_rmse1 = ((yp_pred[:, 0] - y_true[:, 0]) ** 2).mean()
    pool_rmse2 = ((yp_pred[:, 1] - y_true[:, 1]) ** 2).mean()
    pool_rmse = (pool_rmse1 + pool_rmse2) / 2
    pool_rmse = pool_rmse.sqrt().item()

    pool_mape1 = ((yp_pred[:, 0] - y_true[:, 0]) / y_true[:, 0]).abs()
    pool_mape2 = ((yp_pred[:, 1] - y_true[:, 1]) / y_true[:, 1]).abs()
    pool_mape = (pool_mape1 + pool_mape2).sum().item() / 2 / pool_mape1.shape[0] * 100

    print(f'train RMSE: {train_rmse:5.3f}, pool RMSE: {pool_rmse:5.3f}')
    wll = [pool_rmse]
    wmape = [pool_mape]

    end = time.time()
    tTemp.append(end - start)

    # numpy_array = yp_pred.detach().numpy()
    # df = pd.DataFrame(numpy_array)
    # df.to_csv('yp_Pred.csv', index=False)
    # numpy_array = y_pool.numpy()
    # df = pd.DataFrame(numpy_array)
    # df.to_csv('y_pool.csv', index=False)
    # plt.plot(x_train[:, 0].numpy(), y_train[:, 0].numpy(), '.', color='r')
    # plt.plot(x_train[:, 0].numpy(), y_train[:, 1].numpy(), '.', color='r')
    # plt.plot(x_pool[:, 0].numpy(), y_pool[:, 0].numpy(), '-', color='b')
    # plt.plot(x_pool[:, 0].numpy(), y_pool[:, 1].numpy(), '-', color='b')
    # plt.plot(x_pool[:, 0].numpy(), yp_pred[:, 0].detach().numpy(), '--', color='m')
    # plt.plot(x_pool[:, 0].numpy(), yp_pred[:, 1].detach().numpy(), '--', color='m')
    # plt.show()

    # model_filename = f'D:/OneDrive - USTC/Research/Chao/05_Active learning/Code/NN4AL/ModelOfflineTarget{l}.pth'
    # torch.save(ModelTarget.state_dict(), model_filename) # Save model after training with offline target data

    # Online learning
    x_random = torch.tensor(xR.values).float()[:, NStream * l:NStream * l + NStream]#.unsqueeze(1) # load random data
    y_random = torch.tensor(yR.values).float()[:, NStream * l:NStream * l + NStream]

    # Online training
    ModelOnline = ModelTarget
    opt = torch.optim.AdamW(ModelOnline.parameters(), lr=2e-2)

    for epoch in range(NRound-1):
        start = time.time()
        if SelectionMethod==1:
            # Random selection
            new_idxs = [2 * epoch, 2 * epoch +2]
            x_selected = x_random[new_idxs]
            y_selected = y_random[new_idxs]
        else:
            # Active selection
            train_data = TensorFeatureData(x_train)
            pool_data = TensorFeatureData(x_pool)
            new_idxs, _ = select_batch(batch_size=2, models=[ModelOnline],
                                   data={'train': train_data, 'pool': pool_data}, y_train=y_train,
                                   selection_method='lcmd', sel_with_train=True,
                                   base_kernel='grad', kernel_transforms=[('rp', [512])])
            x_selected = x_pool[new_idxs]
            y_selected = y_pool[new_idxs]

        # Merge online data and offline data
        yy = ModelOnline(xx).detach()
        x_train = torch.cat((xx, x_selected), dim=0) # torch.cat((x_train, x_selected), dim=0)
        y_train = torch.cat((yy, y_selected), dim=0) # torch.cat((y_train, y_selected), dim=0)
        # print(y_train)
        # print(y_selected)

        for epoch1 in range(NTOnline):
            y_pred = ModelOnline(x_train)
            loss1 = ((y_pred[:, 0] - y_train[0:, 0]) ** 2).mean()
            loss2 = ((y_pred[:, 1] - y_train[0:, 1]) ** 2).mean()
            loss = (loss1 + loss2) / 2
            train_rmse = loss.sqrt().item()
            loss.backward()
            opt.step()
            opt.zero_grad()

        # NTOnline = NTOnline+50
        # print(NTOnline)
        yp_pred = ModelOnline(x_pool)
        pool_rmse1 = ((yp_pred[:, 0] - y_true[:, 0]) ** 2).mean()
        pool_rmse2 = ((yp_pred[:, 1] - y_true[:, 1]) ** 2).mean()
        pool_rmse = (pool_rmse1 + pool_rmse2) / 2
        pool_rmse = pool_rmse.sqrt().item()

        pool_mape1 = ((yp_pred[:, 0] - y_true[:, 0]) / y_true[:, 0]).abs()
        pool_mape2 = ((yp_pred[:, 1] - y_true[:, 1]) / y_true[:, 1]).abs()
        pool_mape = (pool_mape1 + pool_mape2).sum().item() / 2 / pool_mape1.shape[0] * 100
        print(pool_rmse)
        wll.append(pool_rmse)
        wmape.append(pool_mape)

        end = time.time()
        tTemp.append(end - start)
        # plt.plot(x_pool[:, 0].numpy(), y_pool[:, 0].numpy(), '-', color='b')
        # plt.plot(x_pool[:, 0].numpy(), y_pool[:, 1].numpy(), '-', color='b')
        # plt.plot(x_pool[:, 0].numpy(), yp_pred[:, 0].detach().numpy(), '--', color='m')
        # plt.plot(x_pool[:, 0].numpy(), yp_pred[:, 1].detach().numpy(), '--', color='m')
        # plt.plot(x_train[:, 0].numpy(), y_train[:, 0].numpy(), '.', color='r')
        # plt.plot(x_train[:, 0].numpy(), y_train[:, 1].numpy(), '.', color='r')
        # plt.plot(x_selected[:, 0].numpy(), y_selected[:, 0].numpy(), '.', color='k')
        # plt.plot(x_selected[:, 0].numpy(), y_selected[:, 1].numpy(), '.', color='k')
        # plt.show()

    wl_loss[:, l] = wll
    wl_mape[:, l] = wmape
    ExeTime[:, l] = tTemp

if SelectionMethod==1:
    np.savetxt('Random2D2S2D_RMSE.csv', wl_loss, delimiter=',', header=','.join([f' {i+1}' for i in range(NRep)]), comments='')
    np.savetxt('Random2D2S2D_MAPE.csv', wl_mape, delimiter=',', header=','.join([f' {i + 1}' for i in range(NRep)]),
               comments='')
    np.savetxt('Random2D2S2D_TimeVV.csv', ExeTime, delimiter=',',
               header=','.join([f' {i + 1}' for i in range(NRep)]), comments='')
else:
    np.savetxt('Active2D2S2D_RMSE.csv', wl_loss, delimiter=',', header=','.join([f' {i+1}' for i in range(NRep)]), comments='')
    np.savetxt('Active2D2S2D_MAPE.csv', wl_mape, delimiter=',', header=','.join([f' {i + 1}' for i in range(NRep)]),
               comments='')
    np.savetxt('Active2D2S2D_Time.csv', ExeTime, delimiter=',',
               header=','.join([f' {i + 1}' for i in range(NRep)]), comments='')
