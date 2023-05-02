# %%
import pandas as pd
import pickle
import requests
# %%
with open('./cpddata_dict.pickle','rb') as f:
    cpddata_dict=pickle.load(f)
# %%
# %%
ids=list(cpddata_dict.keys())
df_cpddata=pd.DataFrame([cpddata_dict[id] for id in ids],index=ids,
            columns=['name'])
# %%
# df_cpddata.to_csv('./cpddata.csv')
# %%
data=requests.get('http://rest.kegg.jp/list/cpd')
# %%
data_id=[elem.split('\t')[0] for elem in data.text.splitlines()]
data_name=[elem.split('\t')[1] for elem in data.text.splitlines()]
# %%
data_dict2=dict(zip(data_id,data_name))
# %%
# with open('./cpddata_dict2.pkl','wb') as f:
#     pickle.dump(data_dict2,f)
# %%
data_dict3=dict()
for cpdID in list(data_dict2.keys()):
    data_dict3[cpdID]=data_dict2[cpdID].split(';')[0]
# %%
# with open('./cpddata_dict3.pkl','wb') as f:
#     pickle.dump(data_dict3,f)
# %%
