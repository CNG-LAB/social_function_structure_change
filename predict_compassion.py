'''
This is a python script for running supervised machine learning to predict compassion score

data input : foldname

data output $np.save..... 
'''

import numpy as np
import pandas as pd
import scipy.stats as ss
import sklearn.linear_model as slm
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from mlxtend.feature_selection import SequentialFeatureSelector as SFS
import warnings
warnings.filterwarnings("ignore")

df_data = pd.read_csv ('../data/compassion_data.csv')
sex_code = df_data['sex1'].unique()
if type(sex_code[0]) == str:
    sex_dict = {k: idx for idx, k in enumerate(sex_code)}
    df_data['sex1'] = df_data['sex1'].replace(sex_dict)  


def model_elasticnet(m): #l1_ratio
  print('m', str(m))

  dic = {}
  for i in range(sample): 
    lr = slm.ElasticNetCV(alphas=[0.0001, 0.001, 0.01, 0.1, 1], l1_ratio=m, cv=5)
    sfs = SFS(lr, k_features=7, forward=True, floating=False, n_jobs=-1,
              scoring='neg_mean_absolute_error', cv=False)
    model = sfs.fit(x_correct[i], y_train[i])
    x_1 = model.transform(x_correct[i])
    x_2 = model.transform(x_test_corr[i])
    model2 = lr.fit(x_1, y_train[i])
    y_pred = lr.predict(x_2)
    corr = ss.pearsonr(y_pred, y_test[i])
    a = model.get_metric_dict()[7]
    a['importance'] = model.estimator.coef_
    a['intercept'] = model.estimator.intercept_
    a['alpha_best'] = model.estimator.alpha_
    a['predict_test_r_p'] = np.array(corr)
    a['mean_lr_mse']  = model2.mse_path_.mean()
    dic['train_test_'+str(i+1)] = a
    print('finish model......', str(i+1))
    print(a)
  np.save('../results/'+foldname+'feature_l1ratio_' + str(m)+ '.npy', dic)
  return print('feature_l1ratio_' + str(m)+'   finished')

foldname=['compassion']
for foldname in foldname:

  # IMPORT sample train_test iterations
  sample = 100

  x_train = [None] * sample
  y_train = [None] * sample
  x_test = [None] * sample
  y_test = [None] * sample
  x_correct = [None]  * sample
  x_test_corr = [None] * sample
  for i in range(sample):
    Y_col = 'val1'
    X_cols = df_iris.loc[:, df_iris.columns != Y_col].columns

    x_train[i], x_test[i], y_train[i], y_test[i] = train_test_split(        
          df_data[X_cols], df_data.iloc[:, 1], test_size=0.3, random_state=i)
	
    x_conf = x_train[i].iloc[:,[2,3]]
    y_conf = x_train[i].iloc[:,4:39]
    x, y   = np.array(x_conf), np.array(y_conf)
    model_conf = LinearRegression().fit(x, y)
    y_pred = model_conf.predict(x)
    x_correct[i] = y  - y_pred
    x_conf = x_test[i].iloc[:,[2,3]]
    y_conf = x_test[i].iloc[:,4:39]
    x, y   = np.array(x_conf), np.array(y_conf)
    model_conf = LinearRegression().fit(x, y)
    y_pred = model_conf.predict(x)
    x_test_corr[i] = y  - y_pred

  regulation = [1.0]

  for j in regulation:
    model_elasticnet(m=j)

