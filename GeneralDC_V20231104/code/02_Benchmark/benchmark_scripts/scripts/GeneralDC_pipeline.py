import pandas as pd
import scanpy as sc
import numpy as np
import tangram as tg
import sys


'''
@time 20230823
@desc 1-随机抽取细胞比例; 2-构建模拟spot, 得到pseudo-cell percentage; [3]-将Stage2输入到Stage1。
@update cell-type的数据使用原始的几种类型,不进行log2; spatial的数据使用count
@methods SpatialDC-Train and Auto-supervised model

'''

import sys
import os
import pandas as pd
import scanpy as sc
import numpy as np
from scipy.stats import pearsonr,spearmanr
from sklearn.metrics import mean_squared_error

import matplotlib.pyplot as plt

import torch
import torch.nn as nn
import torch.nn.functional as F

import warnings
warnings.filterwarnings("ignore")
import glob
from sklearn.preprocessing import MinMaxScaler, StandardScaler
import anndata
from scipy.sparse import csr_matrix
from tqdm import tqdm
import math
import random

from datetime import datetime
import torch_geometric
from torch_geometric.nn import VGAE
from torch_geometric.nn import VGAE, GCNConv
from torch_geometric.utils import negative_sampling, remove_self_loops, add_self_loops


random.seed(0)
os.environ['PYTHONHASHSEED'] = '0'
np.random.seed(0)
torch.manual_seed(0)
torch.cuda.manual_seed(0)
torch.cuda.manual_seed_all(0)
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.deterministic = True
torch.use_deterministic_algorithms(True)
os.environ['CUBLAS_WORKSPACE_CONFIG']=':4096:8'

device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')
from scipy.stats import norm
import time
import warnings
warnings.filterwarnings("ignore")



###--------------------------- PART I, Load methods -------------------------------
## 1. Set DNN model
class DNNModel(nn.Module):
    state = {}
    def __init__(self, input_dim, output_dim, drop_rate=0.2):
        super(DNNModel, self).__init__()

        self.inputdim = input_dim
        self.outputdim = output_dim
        self.encoder = nn.Sequential(            
            
            # Layer 1
            nn.Dropout(drop_rate),
            nn.Linear(self.inputdim, self.inputdim // 2),    
            nn.CELU(),

            # Layer 2
            nn.Dropout(drop_rate),            
            nn.Linear(self.inputdim // 2, self.inputdim // 4),# Layer 2            
            nn.CELU(),

            # Layer 3
            nn.Dropout(drop_rate),
            nn.Linear(self.inputdim // 4, self.outputdim),# Layer 3
        )


    def forward(self, x):
        
        z = self.encoder(x)        
        z = F.sigmoid(z)
        z = (z.T / z.sum(axis=1)).T
        

        return z

## 2. Train and Valid 
class PseudoDeconv():        
    def __init__(
        self,        
        num_spots: int=10000,
        test_size=0.1,
        batch_size: int=2048,
        lr: float=0.001,
        lr_decay:float=0.95,
        lr_step=1,
        weight_decay: float=1e-8,        
        start_epoch=1,
        epochs=1,
        adata ="",
        cell_type_obs="celltype",
        
        #simulation data
        sparse_prob:float=0.5,
        maxcell_in_onespot=5,
        print_info = True,
        mean=0,
        sd=5
    ):
        self.cell_type_obs = cell_type_obs
        self.adata = adata

        self.num_spots = num_spots
        self.test_size = test_size    
        self.batch_size = int(batch_size)#int(batch_size * (self.num_spots)/1000)
        self.lr = lr #* np.sqrt(self.batch_size/16)
        self.lr_decay = lr_decay
        self.lr_step = lr_step
        self.weight_decay = weight_decay
        self.start_epoch = start_epoch
        self.epochs = epochs        

        self.celltype_order = self.adata.obs[self.cell_type_obs].cat.categories.values
        # print(self.celltype_order)
        # exit()
        self.maxcell_in_onespot=maxcell_in_onespot

        # self.d_prior = np.ones(self.cell_type)

        self.adata = adata
        self.label_list = self.adata.obs[self.cell_type_obs].cat.categories.tolist()
        self.cell_type = self.adata.obs[self.cell_type_obs].cat.codes.to_numpy()

        self.sparse_prob = sparse_prob
              
        self.print_info = print_info

        self.mean=mean
        self.sd = sd
        self.model = DNNModel(
            input_dim=self.adata.shape[1],            
            output_dim=len(self.label_list),
        ).to(device)

        self.spot_type = ""
        self.spot_data = ""
        # self.model = MyModel(input_size=253, attention_size=253, hidden_size=253, output_size=19).to(device)

    def train(self, retrained=False, pseudo_perc_df = "",std_dev=0.1):
        

        spot_data, spot_type = Utils.generate_simulated_data(sc_data=self.adata,cell_type_obs=self.cell_type_obs, 
                                                                    samplenum=int(self.num_spots))        
                    
        self.spot_type = pd.DataFrame(spot_type, columns=self.label_list)
        self.spot_data = pd.DataFrame(spot_data)
        # print(spot_type)
        train_loader, valid_loader = Utils.prepare_dataloader(spot_data=spot_data, spot_type=spot_type,test_size=self.test_size,batch_size=self.batch_size) # create dataload of sc
        
        model = self.model
        # print(model)
        criterion = torch.nn.MSELoss()#torch.nn.MSELoss()
        # criterion2 = torch.nn.MSELoss()
        params = []
        for key, value in dict(model.named_parameters()).items():
            if 'bias' in key:
                params.append({'params': [value], 'lr': self.lr, 'weight_decay': 0})
            else:
                params.append({'params': [value], 'lr': self.lr, 'weight_decay': self.weight_decay,'betas':[0.9,0.999]})
            
        optimizer = torch.optim.Adam(params)
        
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, self.lr_step, gamma=self.lr_decay)
        # scaler = torch.cuda.amp.GradScaler()
        
        TRAIN_LOSS = []
        VALID_LOSS = []
        VALID_METRIC = []
        # LR = []
        # BEST_VAL = 100       

        for epoch in range(self.start_epoch, self.epochs + 1):

            # ************** Train ************** #
            if(self.print_info):
                print('[Epoch %2d][Lr  %.e]' % (epoch, scheduler.get_last_lr()[0]))
            model.train()
            train_loss = 0

            # Train model.
            for batch_idx, (inputs, targets) in enumerate(train_loader):
                # Send tensors to gpu.
                inputs  = inputs.float().to(device)
                targets = targets.float().to(device)

                # Forward propagation.
                with torch.cuda.amp.autocast():
                    outputs = model(inputs)
                    # outputs = (outputs.T/outputs.sum(axis=1)).T                
                    loss = criterion(targets, outputs)                    
                    
                train_loss += loss.item() #+ loss2.item()

                if math.isnan(train_loss):
                    self.epochs = epoch-1                                        
                    break

                # Backward propagation.
                optimizer.zero_grad()
                loss.backward()
                
                # from GPT3
                # scaler.scale(loss)
                # scaler.unscale_(optimizer)
                optimizer.step() # BP
                # scaler.update() # adjust scale factor
                # scaler.scale(loss).backward()
                # scaler.step(optimizer)
                # scaler.update()

                # if(self.print_info):
                #     print('[Epoch %2d][Iter %4d/ %4d] Loss: %.5f \r'
                #         % (epoch, batch_idx + 1, len(train_loader), train_loss / (batch_idx + 1)), end='')
            # print(targets)
            # print(outputs)
            # exit()
            if math.isnan(train_loss):
                self.epochs = epoch-1
                break
            # Update global variables.
            TRAIN_LOSS.append(train_loss / len(train_loader))            

            # ************** Valid ************** #
            model.eval()
            valid_loss = 0

            # Validate model.
            with torch.no_grad():
                for batch_idx, (inputs, targets) in enumerate(valid_loader):
                    # Send tensors to gpu.
                    inputs  = inputs.float().to(device)
                    targets = targets.float().to(device)

                    # Forward propagation.
                    outputs = model(inputs)             
                    loss = criterion(outputs, targets) #+ criterion(x_recon, inputs)
                    valid_loss += loss.item()
                    
                    # Compute accuracy.
                    # rmse = np.sqrt(mean_squared_error(targets.cpu(), outputs.detach().cpu()))

            # Update global variables.
            VALID_LOSS.append(valid_loss / len(valid_loader))
            # VALID_METRIC.append(rmse)
            if(self.print_info):
                print('[Epoch %2d][Train] Loss: %.5f [Valid] Loss: %.5f ' % (epoch, TRAIN_LOSS[-1], VALID_LOSS[-1]))

            state = {
                'model': model.state_dict(),
                'optimizer': optimizer.state_dict(),
                'epoch': epoch,
                'train_loss': TRAIN_LOSS,
                'valid_loss': VALID_LOSS,
                # 'valid_metric': VALID_METRIC,
                'lr': self.lr,
            }
            model.state = state
            scheduler.step()
            
            # torch.save(state, os.path.join(self.save_dir, f'checkpoint_Brain_CellType_Count_19_Best_'+str(self.epochs)+"_"+str(int(self.num_spots))))

            
        self.model = model
            # torch.optim.swa_utils.update_bn(train_loader, swa_model)
            

    def predict(self, sp_adata):

        inputs = sp_adata.X.toarray() # 2492 * 205

        # print(inputs.shape)
        model = self.model        

        # inputs = np.log1p(inputs)
        inputs = ((inputs.T - inputs.mean(1)) / (inputs.std(1) + 1e-10)).T
        # inputs = (inputs - inputs.mean(0)) / (inputs.std(0) + 1e-10)

        # ss = StandardScaler()
        # inputs = ss.fit_transform(inputs.T).T
        # ss_test_x = ss.fit_transform(test_x.T).T

        # print(inputs)
        # exit()
        inputs = torch.from_numpy(inputs)

        model.eval()
        inputs  = inputs.float().to(device)
        outputs= model(inputs)#
        outputs = outputs.detach().cpu().numpy()               

        preds = pd.DataFrame(outputs)
        preds.columns = self.celltype_order 

        return preds
    


## Utils
class Utils():

    def __init__(self) -> None:
        pass
    # for stage 1
    def prepare_dataloader(spot_data, spot_type,test_size,batch_size):
        # Set training and validation sets.
        num_spots = spot_data.shape[0]
        if test_size < 1:
            test_size = int(num_spots * test_size)
        _train_set = scale_dataset(spot_data[test_size:], spot_type[test_size:])
        _valid_set = scale_dataset(spot_data[:test_size], spot_type[:test_size])
        # print('  Train: %d spots,  Valid: %d spots' % (len(_train_set), len(_valid_set)))

        _train_loader = torch.utils.data.DataLoader(_train_set, batch_size=batch_size, shuffle=True, pin_memory=True)
        _valid_loader = torch.utils.data.DataLoader(_valid_set, batch_size=batch_size * 2, shuffle=False, pin_memory=True)
        return _train_loader, _valid_loader
         
    def generate_simulated_data(sc_data, cell_type_obs="celltype",samplenum=100):

        # print(sc_data)
        if isinstance(sc_data.X, np.ndarray):
            pass
        else:
            sc_data.X = sc_data.X.toarray() #73 * 10466


        num_celltype = len(sc_data.obs[cell_type_obs].value_counts())
                    
        prop = np.random.dirichlet(np.ones(num_celltype), samplenum)

        spot_data = np.matmul(prop,sc_data.X )

        return (spot_data, prop)
    
class scale_dataset(torch.utils.data.Dataset):

    def __init__(self, inputs, labels):
        super(scale_dataset, self).__init__()
        # self.inputs = inputs
        # inputs = np.log1p(inputs) # log2 normalization        
        self.inputs = ((inputs.T - inputs.mean(1)) / (inputs.std(1) + 1e-10)).T
        
        self.labels = labels

    def __getitem__(self, index):
        _input = self.inputs[index]
        _label = self.labels[index]        
        return _input, _label

    def __len__(self):
        return len(self.labels)    




# run stage 1 and 2
def run_GeneralDC(sc_adata,sp_adata,output_file_path, sc_epochs,weight_decay=1e-4,print_info=False,batch_size=128,lr=0.001,num_spot=10000):
    #Stage 1
    
    intersect = np.intersect1d(sc_adata.var_names, sp_adata.var_names)
    sc_adata = sc_adata[:, intersect].copy()
    sp_adata = sp_adata[:, intersect].copy()

    obs_index = sp_adata.obs.index.values

    scDeconv_model = PseudoDeconv(adata=sc_adata,epochs=sc_epochs,print_info=print_info,num_spots=num_spot, lr=lr,
                                weight_decay= weight_decay,batch_size=batch_size)
    scDeconv_model.train()
    #plot sc
    y = scDeconv_model.model.state["train_loss"]  # 定义y轴为从0到1连续取值的数组
    valid_loss = scDeconv_model.model.state["valid_loss"]
    plt.plot(np.arange(0,scDeconv_model.epochs), y)  # 绘制曲线图
    plt.xlabel('Epochs')
    plt.ylabel('Train Loss')
    plt.title(f'Simu_Spot:[{num_spot}] | Train Loss:[{"%.4f" % y[-1]}] | Valid Loss:[{"%.4f" % valid_loss[-1]}]')
    plt.savefig(f"{output_file_path}.png")


    # checkpoint = torch.load(model_path, map_location=torch.device('cpu'))
    # scDeconv_model.model.load_state_dict(checkpoint['model'])
    # scDeconv_model.model.eval()  

    psedu_pred_df = scDeconv_model.predict(sp_adata) # use intersect genes of spatial_adata

    psedu_pred_df.index = obs_index
    psedu_pred_df.to_csv(f'{output_file_path}.csv')

    return(psedu_pred_df)


###--------------------------- PART II, Set Parms ----------------------------------
sc_epoch = 50
num_spot = 50000

sc_file_path = sys.argv[1]
spatial_file_path = sys.argv[2]
celltype_key = sys.argv[3]
output_file_path = sys.argv[4]

sc_adata = sc.read_h5ad(sc_file_path)
sp_adata = sc.read_h5ad(spatial_file_path)

###--------------------------- PART III, Train model and predict--------------------
# Stage 1,2, Train and Predict 
run_GeneralDC(sc_adata=sc_adata, sp_adata=sp_adata, sc_epochs=sc_epoch,num_spot=num_spot,print_info=True,
                output_file_path=output_file_path)
