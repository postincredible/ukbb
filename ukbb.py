

import os
import pandas as pd
import numpy as np

pth=os.getcwd()
rel_var_path=pth.split('code')[0]+'disease/'
rel_var_path

def load_data_by_fid(fid):
    df_tab1_i0_comp=pd.read_csv('/oasis/scratch/comet/yhuan162/temp_project/ukbb/data/i0/ukb22598_i0_comp.csv')

    if int(fid) in df_tab1_i0_comp.fid.values.tolist():
        fid_num=fid
        
        var_type=df_tab1_i0_comp[df_tab1_i0_comp['fid']==int(fid_num)].Type.values[0]

        var_type_list=['con','cur','dat','int','tex','tim','cas','cam']
        var_type_list_full=['Continuous','Curve','Date','Integer','Text','Time','Categorical (single)', 'Categorical (multiple)']

        path_p1='/oasis/scratch/comet/yhuan162/temp_project/ukbb/data/i0/var_'

        if var_type in var_type_list_full:
            vtyp=var_type_list[var_type_list_full.index(var_type)]

        loadpath=path_p1+str(vtyp)+'/'
        os.chdir(path_p1+str(vtyp))
        list_folder=os.listdir() 

        pname1=str(vtyp)+str(fid_num)+'i0.csv'
        pname2='vec_'+str(vtyp)+str(fid_num)+'i0.csv'

        if pname1 in list_folder:
            print('fid ' + str(fid_num) + ' is a single-measure '+str(var_type).lower()+' variable')
            fpname=list_folder[list_folder.index(pname1)]
            df_load=pd.read_csv(loadpath+fpname)

        elif pname2 in list_folder:
            print('fid ' + str(fid_num) + ' is a multiple-measure '+str(var_type).lower()+' variable')
            fpname=list_folder[list_folder.index(pname2)]
            df_load=pd.read_csv(loadpath+fpname, sep='\t')
        return df_load
    
    else:
        print('fid not found, please try again')

        
df_tab1_i0_comp=pd.read_csv('/oasis/scratch/comet/yhuan162/temp_project/ukbb/data/i0/ukb22598_i0_comp.csv')

def cd_path(path):
    """
    Check if path exists, if not create it
    """
    if os.path.exists(path):
        print(path+' already exists!')
    else:
        os.mkdir(path)
        print(path+' is now created!')

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
            
    savepath1=rel_var_path                        ##### CHANGE path if needed
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










################ functions in std UKBB data ######################

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
    
    


def gen_event_ind_from_list(fid_int, event_code_list):
    """
    return a dataframe that contains indicator variables for each pair of event and its related ICD code

    """
    df_pool=[]
    for lev1 in event_code_list:
        print('load event: '+ str(lev1[0]))
        for_event= lev1[0]
        for_code= lev1[1]
        #df_name= 'df_ind'+str(fid_int)+'_'+str(for_event)
        df_pool_element=mm_gen_ind_list(fid_in=fid_int, evt=for_event, key_code_list=for_code)
        df_pool.append(df_pool_element)
    df_pooled=pd.concat(df_pool,axis=1)
    df_pooled=df_pooled.loc[:,~df_pooled.columns.duplicated()]  # drop duplicated 'eid' columns
    return df_pooled



def gen_event_ind_from_multi_var(fid_list, event_code_list,detail=False):
    """
    return a dataframe that contains indicator variables for each event combined multiple icd measurements

    """
    f_pool=[]
    for f in fid_list:
        print('\n working on fid= '+str(f))
        f_pool_element= gen_event_ind_from_list(fid_int=f, event_code_list=event_code_list)
        f_pool.append(f_pool_element)
    f_pooled= pd.concat(f_pool,axis=1)
    f_pooled=f_pooled.loc[:,~f_pooled.columns.duplicated()]  # drop duplicated 'eid' columns
    
    if detail==True:
        return f_pooled
    
    if detail==False:
        ind_pool=[]
        for e in event_code_list:
            event=e[0]
            df_pre=f_pooled.filter(regex=event)
            ind_name='icd_ind_'+str(event)
            leid=f_pooled.eid
            ind_sum=df_pre.sum(axis=1)
            df_e=pd.DataFrame({'eid':leid, ind_name:ind_sum})
            df_e[ind_name]=df_e[ind_name].apply(lambda y: 1 if y>0 else y)
            df_e=df_e.loc[:,~df_e.columns.duplicated()]  # drop duplicated 'eid' columns
            ind_pool.append(df_e)
            
            #df_e=f_pooled[['eid']].copy()
            #ind_name='icd_ind_'+str(event)
            #df_e[str(ind_name)]=df_pre.sum(axis=1)
            #df_e.ind_name=df_pre.sum(axis=1)
            #df_e[ind_name]=df_pre.sum(axis=1)

            #df_e.ind_name=df_e.ind_name.apply(lambda y: 1 if y>0 else y)
            #ind_pool.append(df_e)
        ind_pooled= pd.concat(ind_pool,axis=1)
        ind_pooled=ind_pooled.loc[:,~ind_pooled.columns.duplicated()]  # drop duplicated 'eid' columns
        return ind_pooled
    
    

    

############# functions for HES ################

def hes_gen_ind_raw(icd,hesin_dfin,key_code,evnt, detail=False):
    """
    return a dataframe that contains indicator variable for a specific 'key_code' in HES data
        use 'detail= True' to get the detail matched code info
    """

    #dfc=load_data_by_fid(fid_int)
    #df_icd9m=dfc.copy()
    #dfa=hesin[['eid',str(icd)]].copy()
    dfa=hesin_dfin[['eid','record_id',str(icd)]].copy()


    
    dfa_lst=dfa[dfa.columns[dfa.columns.get_loc(str(icd))]].values.tolist()
    pre0=lst_ind(dfa_lst,str(key_code))

    gen_hes_name='hes_'+str(icd)+'_'+str(evnt)+str(key_code)
    gen_ind_name='ind_'+str(icd)+'_'+str(evnt)+str(key_code)

    dfa[str(gen_hes_name)]=pre0
    dfa[dfa.columns[dfa.columns.get_loc(str(gen_hes_name))]]=dfa[dfa.columns[dfa.columns.get_loc(str(gen_hes_name))]].apply(lambda y: np.nan if len(y)==0 else y )
    
    dfa[str(gen_ind_name)]=pre0
    dfa[dfa.columns[dfa.columns.get_loc(str(gen_ind_name))]]=dfa[dfa.columns[dfa.columns.get_loc(str(gen_ind_name))]].apply(lambda y: 0 if len(y)==0 else 1 )
    
    print('\nHES '+str(icd)+' ',str(evnt)+'('+str(key_code)+')'+' count: '+str(dfa[dfa.columns[dfa.columns.get_loc(str(gen_hes_name))]].count())+',\nFreq_tab \n'+str(dfa[dfa.columns[dfa.columns.get_loc(str(gen_ind_name))]].value_counts()))
    dfb=dfa[['eid','record_id',str(gen_ind_name)]]   
    
    if detail==True:
        return dfa
    else:
        return dfb




def hes_gen_ind_list(icd_in, hesin_dfin, key_code_list, evt, detai=False):
    """
    return a dataframe that contains indicator variables for each specific 'key_code' in 'key_code_list'
        use 'detai= True' to get the detail matched codes info
    """
    dfcl=[]
    if detai==False:
        for l in key_code_list:
            df_l=hes_gen_ind_raw(icd_in,hesin_dfin, l, str(evt), detail=False)
            dfcl.append(df_l)
        dfcl_merge=pd.concat(dfcl,axis=1)
        dfcl_merge=dfcl_merge.loc[:,~dfcl_merge.columns.duplicated()]  # drop duplicated 'eid' columns
        return dfcl_merge

    if detai==True:
        for l in key_code_list:
            df_l=hes_gen_ind_raw(icd_in,hesin_dfin, l, str(evt), detail=True)
            dfcl.append(df_l)
        dfcl_merge=pd.concat(dfcl,axis=1)
        dfcl_merge=dfcl_merge.loc[:,~dfcl_merge.columns.duplicated()]  # drop duplicated 'eid' columns
        return dfcl_merge

    

def hes_gen_event_ind_from_list(icd_var, hes_df, event_code_list):
    """
    return a dataframe that contains indicator variables for each pair of event and its related ICD code from HES database

    """
    df_pool=[]
    for lev1 in event_code_list:
        print('load event: '+ str(lev1[0]))
        for_event= lev1[0]
        for_code= lev1[1]
        #df_name= 'df_ind'+str(fid_int)+'_'+str(for_event)
        df_pool_element=hes_gen_ind_list(icd_in=icd_var,hesin_dfin=hes_df, evt=for_event, key_code_list=for_code)
        df_pool.append(df_pool_element)
    df_pooled=pd.concat(df_pool,axis=1)
    df_pooled=df_pooled.loc[:,~df_pooled.columns.duplicated()]  # drop duplicated 'eid' columns
    return df_pooled


def hes_gen_event_ind_from_multi_var(icd_var_list, hes_dfin, event_code_list,detail=False):
    """
    return a dataframe that contains indicator variables for each event combined multiple icd measurements

    """
    f_pool=[]
    for f in icd_var_list:
        print('\n working on icd_var= '+str(f))
        f_pool_element= hes_gen_event_ind_from_list(icd_var=f, hes_df=hes_dfin, event_code_list=event_code_list)
        f_pool.append(f_pool_element)
    f_pooled= pd.concat(f_pool,axis=1)
    f_pooled=f_pooled.loc[:,~f_pooled.columns.duplicated()]  # drop duplicated 'eid' columns
    
    if detail==True:
        return f_pooled
    
    if detail==False:
        ind_pool=[]
        for e in event_code_list:
            event=e[0]
            df_pre=f_pooled.filter(regex=event)
            ind_name='hes_icd_ind_'+str(event)
            leid=f_pooled.eid
            lrec=f_pooled.record_id
            ind_sum=df_pre.sum(axis=1)
            df_e=pd.DataFrame({'eid':leid,'record_id':lrec,ind_name:ind_sum})
            df_e[ind_name]=df_e[ind_name].apply(lambda y: 1 if y>0 else y)
            df_e=df_e.loc[:,~df_e.columns.duplicated()]  # drop duplicated 'eid' columns
            ind_pool.append(df_e)
        ind_pooled= pd.concat(ind_pool,axis=1)
        ind_pooled=ind_pooled.loc[:,~ind_pooled.columns.duplicated()]  # drop duplicated 'eid' columns
        return ind_pooled         
            

