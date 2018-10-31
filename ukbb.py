import os
import pandas as pd
import numpy as np


os.chdir('/projects/ps-janssen3/dsci-pa/yhuan162/temp_project/all_codes')
from ukbb_ldbf import load_data_by_fid

df_tab1_i0_comp=pd.read_csv('/projects/ps-janssen3/dsci-pa/yhuan162/temp_project/ukbb/data/i0/ukb22598_i0_comp.csv')

def chk_unique_eid(df):
    """
    Check unique eid number for a dataframe
    """
    print('loaded df has unique eid count: '+ str(len(df.eid.unique())))

    
    
def search_des(keyword):
    """
    search 'keyword' related variable based on the variable description
    """
    klow=str(keyword).lower()
    df_tab1_i0_comp['des']=df_tab1_i0_comp.Description.str.lower()
    key_des=df_tab1_i0_comp[df_tab1_i0_comp.des.str.contains(klow)][['fid','Type','obs_ct','Description','DC']]
    return key_des


def related_vars(key_list, dis):
    """
    return a dataframe contains all searched 'keyword' related variable in 'key_list'
    """
    savepath1='/oasis/scratch/comet/yhuan162/temp_project/ukbb/data/disease/'
    savepath2=savepath1+str(dis).upper()
    
    if os.path.exists(savepath2):
        os.chdir(savepath2)
        d_lst=[]
        for k in key_list:
            df_k=search_des(str(k).strip())
            d_lst.append(df_k)

        d_coma=pd.concat(d_lst)
        d_comb=d_coma.drop_duplicates()
        print('Searched keyword(s): '+str(key_list)+'\n'+'save '+str(dis)+'_related_vars_chk.csv file at '+str(savepath2))
        filename=str(dis)+'_related_vars_chk.csv'
        d_comb.to_csv(filename, index=None)
        return d_comb
    
    else: 
        os.mkdir(savepath2)
        os.chdir(savepath2)
        d_lst=[]
        for k in key_list:
            df_k=search_des(str(k).strip())
            d_lst.append(df_k)

        d_coma=pd.concat(d_lst)
        d_comb=d_coma.drop_duplicates()
        print('Searched keyword(s): '+str(key_list)+'\n'+'save '+str(dis)+'_related_vars_chk.csv file at '+str(savepath2))
        filename=str(dis)+'_related_vars_chk.csv'
        d_comb.to_csv(filename, index=None)
        return d_comb




def lst_ind(dfa_list,ind_val):
    """
    return a list of icd code that match with 'ind_val'
    """
    pre0=[]
    for i in dfa_list:
        if pd.isnull(i):
            pre0.append([])
        elif pd.notnull(i):
            si=[]
            jl=i.split(',')
            for ei in jl:
                ef=ei.replace(',','')
                efa,efb,efc=ef.partition(str(ind_val))
                if efa=='':
                    si.append(ef)
            pre0.append(si)
    return pre0


def mm_gen_ind_raw(fid_int,key_code,evnt, detail=False, get_ct=False, ct_only=False):
    """
    return a dataframe that contains indicator variable for a specific 'key_code' in UKBB std data
        use 'detail=True' to get the detail matched code info
        use 'get_ct=True' to get the count for matched code
        use 'ct_only=True' to return count only
    """

    dfc=load_data_by_fid(fid_int)
    #df_icd9m=dfc.copy()
    dfa=dfc.copy()

    dfa_lst=dfa[dfa.columns[1]].values.tolist()
    
    pre0=lst_ind(dfa_lst,str(key_code))
    
    gen_fid_name='fid'+str(fid_int)+'_'+str(evnt)+str(key_code)
    gen_ind_name='ind'+str(fid_int)+'_'+str(evnt)+str(key_code)
    gen_count_name='count'+str(fid_int)+'_'+str(evnt)+str(key_code)
    
    dfa[str(gen_fid_name)]=pre0
    dfa[dfa.columns[dfa.columns.get_loc(str(gen_fid_name))]]=dfa[dfa.columns[dfa.columns.get_loc(str(gen_fid_name))]].apply(lambda y: np.nan if len(y)==0 else y )
    
    dfa[str(gen_ind_name)]=pre0
    dfa[dfa.columns[dfa.columns.get_loc(str(gen_ind_name))]]=dfa[dfa.columns[dfa.columns.get_loc(str(gen_ind_name))]].apply(lambda y: 0 if len(y)==0 else 1 )
    
    dfa[str(gen_count_name)]=pre0
    dfa[dfa.columns[dfa.columns.get_loc(str(gen_count_name))]]=dfa[dfa.columns[dfa.columns.get_loc(str(gen_count_name))]].apply(lambda y: 0 if len(y)==0 else len(y) )
    
    print('fid '+str(fid_int)+' ',str(evnt)+str(key_code)+' count: '+str(dfa[dfa.columns[dfa.columns.get_loc(str(gen_fid_name))]].count())+' ind from '+str(dfa[dfa.columns[dfa.columns.get_loc(str(gen_ind_name))]].count()))
    dfb=dfa[['eid',str(gen_ind_name),str(gen_count_name)]]
    #dfb=dfa[['eid',str(gen_ind_name)]]
    
    if ct_only==False:
        if detail==True:
            if get_ct==True:
                return dfa
            if get_ct==False:
                return dfa.drop([str(gen_count_name)],axis=1)
        else:
            if get_ct==True:
                return dfb
            if get_ct==False:
                return dfb.drop([str(gen_count_name)],axis=1)
        
    if ct_only==True:
        return dfb.drop([str(gen_ind_name)],axis=1)

    
    
        
def mm_gen_ind_list(fid_in, key_code_list, evt, detai=False, get_ct=False, ct_only=False):
    """
    return a dataframe that contains indicator variables for each specific 'key_code' in 'key_code_list'
        use 'detai= True' to get the detail matched codes info
        use 'get_ct=True' to get the count for matched codes
        use 'ct_only=True' to return counts only
    """
    dfcl=[]
    
    if ct_only==False:
    
        if detai==False:
            if get_ct==False:
                for l in key_code_list:
                    df_l=mm_gen_ind_raw(fid_in, l, str(evt), detail=False, get_ct=False, ct_only=False)
                    dfcl.append(df_l)
                dfcl_merge=pd.concat(dfcl,axis=1)
                dfcl_merge=dfcl_merge.loc[:,~dfcl_merge.columns.duplicated()]  # drop duplicated 'eid' columns
                return dfcl_merge
        
            if get_ct==True:
                for l in key_code_list:
                    df_l=mm_gen_ind_raw(fid_in, l, str(evt), detail=False, get_ct=True, ct_only=False)
                    dfcl.append(df_l)
                dfcl_merge=pd.concat(dfcl,axis=1)
                dfcl_merge=dfcl_merge.loc[:,~dfcl_merge.columns.duplicated()]  # drop duplicated 'eid' columns
                return dfcl_merge
        
        
        
        if detai==True:
            if get_ct==False:
                for l in key_code_list:
                    df_l=mm_gen_ind_raw(fid_in, l, str(evt), detail=True, get_ct=False, ct_only=False)
                    dfcl.append(df_l)
                dfcl_merge=pd.concat(dfcl,axis=1)
                dfcl_merge=dfcl_merge.loc[:,~dfcl_merge.columns.duplicated()]  # drop duplicated 'eid' columns
                return dfcl_merge
        
            if get_ct==True:
                for l in key_code_list:
                    df_l=mm_gen_ind_raw(fid_in, l, str(evt), detail=True, get_ct=True, ct_only=False)
                    dfcl.append(df_l)
                dfcl_merge=pd.concat(dfcl,axis=1)
                dfcl_merge=dfcl_merge.loc[:,~dfcl_merge.columns.duplicated()]  # drop duplicated 'eid' columns
                return dfcl_merge

    if ct_only==True:
        for l in key_code_list:
            df_l=mm_gen_ind_raw(fid_in, l, str(evt), detail=False, get_ct=False, ct_only=True)
            dfcl.append(df_l)
        dfcl_merge=pd.concat(dfcl,axis=1)
        dfcl_merge=dfcl_merge.loc[:,~dfcl_merge.columns.duplicated()]  # drop duplicated 'eid' columns
        return dfcl_merge




def cal_vec_con(df, get_list):
    """
    return a dataframe that contains processed results from the 'get_list' request for variable in vector(numberic) format
         - Currently, output results can be 'min','max','median','mean'.
    """

    cal_list = ['min','max','median','mean']
    func_list=[np.min, np.max, np.median, np.mean]
    var_list=['vec_min','vec_max','vec_median','vec_mean']
    
    for e in get_list:
        if e in cal_list:
            get_type = cal_list[cal_list.index(e)]
            l = df.iloc[:,1].values.tolist()
            ln = []
            for i in l:
                if pd.isnull(i):
                    ln.append(i)
                elif pd.isnull(i)==False:
                    jl=[float(j) for j in i.split(',')]
                    jm=func_list[cal_list.index(e)](jl)
                    ln.append(jm)
            df[var_list[cal_list.index(e)]]=ln
        if e not in cal_list:
            print(str(e)+' not calculated')
    return df