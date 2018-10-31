import os
import pandas as pd
import numpy as np

def load_data_by_fid(fid):
    '''
    return a dataframe that has the eid and the 'fid' variable  
    '''
    df_tab1_i0_comp=pd.read_csv('/projects/ps-janssen3/dsci-pa/yhuan162/temp_project/ukbb/data/i0/ukb22598_i0_comp.csv')

    if int(fid) in df_tab1_i0_comp.fid.values.tolist():
        fid_num=fid
        
        var_description = df_tab1_i0_comp[df_tab1_i0_comp['fid']==int(fid_num)].Description.values[0]
        var_type=df_tab1_i0_comp[df_tab1_i0_comp['fid']==int(fid_num)].Type.values[0]

        var_type_list=['con','cur','dat','int','tex','tim','cas','cam']
        var_type_list_full=['Continuous','Curve','Date','Integer','Text','Time','Categorical (single)', 'Categorical (multiple)']

        path_p1='/projects/ps-janssen3/dsci-pa/yhuan162/temp_project/ukbb/data/i0/var_'

        if var_type in var_type_list_full:
            vtyp=var_type_list[var_type_list_full.index(var_type)]

        loadpath=path_p1+str(vtyp)+'/'
        os.chdir(path_p1+str(vtyp))
        list_folder=os.listdir() 

        pname1=str(vtyp)+str(fid_num)+'i0.csv'
        pname2='vec_'+str(vtyp)+str(fid_num)+'i0.csv'

        if pname1 in list_folder:

            print('fid ' + str(fid_num) + ' is a single-measure '+str(var_type).lower()+' variable, which is \n'+str(var_description))
            fpname=list_folder[list_folder.index(pname1)]
            df_load=pd.read_csv(loadpath+fpname)

        elif pname2 in list_folder:

            print('fid ' + str(fid_num) + ' is a single-measure '+str(var_type).lower()+' variable, which is \n'+str(var_description))
            fpname=list_folder[list_folder.index(pname2)]
            df_load=pd.read_csv(loadpath+fpname, sep='\t')
        return df_load
    
    else:
        print('fid not found, please try again')
