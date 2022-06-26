# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 21:17:45 2022

@author: keval
~~~~ Final version of code for extracting Drugs from Clinical Trials
"""

print("Starting the program")

import csv     # imports the csv module
import sys      # imports the sys module
import re
import os
import pathlib
import numpy as np
import pandas as pd

print("About to load the files")

filesIssue = False

curPath = os.path.abspath(".")
print("curPath: ",curPath)

# We make the assumption that this code is called from the repositorytechnical_test_clean
# alt_path_clinical_trials_2015 = curPath+"\\technical_test_clean\clinical_trials_2015.jsonl"
# alt_path_drugs = curPath+"\\technical_test_clean\drugs.csv"
# alt_path_usan_stems_namefix = curPath+"\\technical_test_clean\\usan_stems_namefix.csv"
alt_path_clinical_trials_2015 = curPath+"\clinical_trials_2015.jsonl"
alt_path_drugs = curPath+"\drugs.csv"
alt_path_usan_stems_namefix = curPath+"\\usan_stems_namefix.csv"

print("alt_path_clinical_trials_2015 : ",alt_path_clinical_trials_2015 )
print("alt_path_drugs : ",alt_path_drugs )
print("alt_path_usan_stems_namefix : ",alt_path_usan_stems_namefix )

# # # -- Load the files with pandas
try:
    f_trials_2015 = pd.read_json('data/clinical_trials_2015.jsonl', lines=True, encoding="utf8") # opens the jsonl file # 
    f_drugs = pd.read_csv('data/drugs.csv', encoding="utf8") # opens the csv file # 
    f_usan_stems = pd.read_csv('data/usan_stems_namefix.csv', encoding="utf8") # opens the csv file # 
except ValueError as e:
    print("Error trying to open a file. ValueError: ",e)
    try:
        f_trials_2015 = pd.read_json(alt_path_clinical_trials_2015, lines=True, encoding="utf8") # opens the jsonl file # 
        f_drugs = pd.read_csv(alt_path_drugs, encoding="utf8") # opens the csv file # 
        f_usan_stems = pd.read_csv(alt_path_usan_stems_namefix, encoding="utf8") # opens the csv file # 
    except ValueError as e:
        print("Error trying to open a file. ValueError: ",e)
        filesIssue = True

if filesIssue: 
    quit()

print("loaded the files")

if not os.path.exists('data'):
    os.makedirs('data')

print("Set a directory data to save the outputs")


# ~~~~~~~~ for task 1
print("Task 1")


# ---- Change drugs.csv into a dictionary from altLabel_list into its itemLabel. 
def prepare_dict_drugs(f_drugs):
    dict_drugs = {}
    for index,r in f_drugs.iterrows():
        cur_item = f_drugs.iloc[index,:].itemLabel
        cur_alt_item_list = f_drugs.iloc[index,:].altLabel_list.split('|')
        # ~~ We make the choice to consider that the name of the medicine can be set in the list of alternative names in all cases
        if not(cur_item in cur_alt_item_list): cur_alt_item_list.append(cur_item)
        len_alt_list_item = len(cur_alt_item_list)
        np_cur_alt_item_list = np.zeros(len_alt_list_item, dtype=np.dtype('U25'))
        np_cur_alt_item_list.fill(cur_item)
        cur_item_arr = np_cur_alt_item_list
        d_tst = dict(zip(cur_alt_item_list,cur_item_arr))
        dict_drugs.update(d_tst)
    return dict_drugs

dict_drugs = prepare_dict_drugs(f_drugs)
# we also need a dict without spaces
dict_drugs_no_space = {str(k).strip(): str(v).strip() for k, v in dict_drugs.items()}


# Prepare the output for task 1
def set_struct_task1 (f_trials_2015):
    l_task1 = []
    f_trials_2015 = f_trials_2015.sort_values("nct_id")
    for index,row in f_trials_2015.iterrows():
        if (f_trials_2015.iloc[index,:].intervention_type == "Drug"):
            cur_drugs_arr = f_trials_2015.iloc[index,:].intervention_name.split('+')
            obj_nct_id = {"nct_id":f_trials_2015.iloc[index,:].nct_id} #, "drugs": cur_drugs_arr}            
            matching_drugs = []
            for i in cur_drugs_arr:
                i = i.strip()
                if(i in dict_drugs_no_space):
                    matching_drugs.append(i)    
            rem_space = list(map(lambda x : x.strip(),cur_drugs_arr))
            obj_nct_id["drugs"] = matching_drugs    
            empty_arr = len(l_task1)==0        
            if len(matching_drugs)>0:
                if empty_arr:
                        l_task1.append(obj_nct_id)
                else:
                    match_nct_id = l_task1[len(l_task1)-1]["nct_id"] == obj_nct_id["nct_id"]
                    if match_nct_id:
                        l_task1[len(l_task1)-1]["drugs"].append(obj_nct_id["drugs"][0]) # need to append to existing nct_id.
                    else:
                        l_task1.append(obj_nct_id)                        
    return l_task1    

l_task1 = set_struct_task1(f_trials_2015)

# ~~ write the result. Seems fine now?
pd_task_output_1 = pd.DataFrame(l_task1)
#create JSON file 
json_file = pd_task_output_1.to_json(orient='records') 
#export JSON file
with open('data/pd_task_output_1.json', 'w') as f:
    f.write(json_file)

# ~~~~ for task 2
# for each drug present in pd_task_output_1, we want to label definitions of f_usan_stems
print("Task 2")
# Returns True if all elements of list are None
def verify_list_allNone(l):
    l_filt = list ( filter( None, l) )
    return len(l_filt)==0

# Dirty approach with loops, but we are in a rush
def verify_matches_struct(l_task1 ,usans_stems):
    result = []
    cntrT = 0
    panda_t1 = pd.DataFrame(l_task1)
    alt_panda_t1 = panda_t1.explode("drugs")
    all_drugs = list(set(list(alt_panda_t1["drugs"])))
    all_drugs.sort()
    all_drugs = ' '.join(all_drugs).split()
        
    cntrLoop = 0
    for drug in all_drugs:
        # print("====== drug: ",drug,", cntrLoop: ",cntrLoop)
        if cntrLoop < len(all_drugs):
            cntrLoop+=1
            for index, row in f_usan_stems.iterrows():
                curName = row["name"]
                curStem = row["stem"]
                isNanName = False
                isStr = type(row["name"]) == str
                if not isStr:
                    isNanName = np.isnan(row["name"])
                isSubGroup = False
                if not isNanName:
                    isSubGroup = "subgroup" in row["name"]
                isSubClass = isNanName or isSubGroup
                curClassType = "class"
                if isSubClass:
                    curClassType="subClass"

                regexPattern = ""
                if isinstance(curStem, str):
                    curStem = curStem.replace("=","")
                    stemRegex = curStem.replace("-","")
                    checkAny=False
                    checkBeg=False
                    checkEnd=False
                    if curStem[0] == "-" and curStem[len(curStem)-1] == "-": # Any
                        checkAny=True
                        regexPattern = ".*"+stemRegex+".*" 
                    elif curStem[0] == "-" == "-": # End
                        checkEnd=True
                        regexPattern = ".*"+stemRegex+"$"
                    elif curStem[len(curStem)-1] == "-": # Beg
                        checkBeg=True
                        regexPattern = "^"+stemRegex
                    
                    drug_match = drug.replace(" ","") 
                    
                    if (checkAny or checkBeg or checkEnd) and (re.search( regexPattern , drug_match)):
                        if len(result) == 0:
                            result.append({
                                "drug":drug,
                                "usan_codes":[
                                    {
                                        "description": row["definition"],
                                        "type":curClassType
                                    }                                    
                                ]
                            })
                        else:
                            if result[len(result)-1]["drug"] == drug:
                                # verify that the description isn't here already.
                                definition_already_in = any(pd.DataFrame(result[len(result)-1]["usan_codes"])["description"].str.contains(row["definition"],regex=False))
                                if not definition_already_in:
                                    result[len(result)-1]["usan_codes"].append({
                                        "description": row["definition"],
                                        "type":curClassType                            
                                    })
                            else:
                                result.append({
                                    "drug":drug,
                                    "usan_codes":[
                                        {
                                            "description": row["definition"],
                                            "type":curClassType
                                        }
                                    ]
                                })
        else:
            break
    return result

pd_task2 = verify_matches_struct(l_task1, f_usan_stems)

# write the output pd_task2
# ~~ write the result. Seems fine now?
pd_task_output_2 = pd.DataFrame(pd_task2)
#create JSON file 
json_file2 = pd_task_output_2.to_json(orient='records') 
#export JSON file
with open('data/pd_task_output_2.json', 'w') as f:
    f.write(json_file2)

# ~~~~ For task 3: Generate counts of trials by USAN class
print("Task 3")
# We aim to make a similar approach than previously but with differences to list trials according to that

# alt_panda_t1
panda_pd_task2 = pd.DataFrame(pd_task2)
alt_panda_pd_task2 = panda_pd_task2.explode("usan_codes")


def set_drugs_classes(l_task1 ,usans_stems):
    result = []
    cntrT = 0
    panda_t1 = pd.DataFrame(l_task1)
    alt_panda_t1 = panda_t1.explode("drugs")
    all_drugs = list(set(list(alt_panda_t1["drugs"])))
    all_drugs.sort()
    all_drugs = ' '.join(all_drugs).split()
        
    cntrLoop = 0
    for drug in all_drugs:
        # print("======")
        # print("drug: ",drug,", cntrLoop: ",cntrLoop)
        if cntrLoop < len(all_drugs):
            cntrLoop+=1
            for index, row in f_usan_stems.iterrows():
                curName = row["name"]
                curStem = row["stem"]
                isNanName = False
                isStr = type(row["name"]) == str
                if not isStr:
                    isNanName = np.isnan(row["name"])
                isSubGroup = False
                if not isNanName:
                    isSubGroup = "subgroup" in row["name"]
                isSubClass = isNanName or isSubGroup
                curClassType = "class"
                if isSubClass:
                    curClassType="subClass"

                regexPattern = ""
                if isinstance(curStem, str):
                    curStem = curStem.replace("=","")
                    stemRegex = curStem.replace("-","")
                    checkAny=False
                    checkBeg=False
                    checkEnd=False
                    if curStem[0] == "-" and curStem[len(curStem)-1] == "-": # Any
                        checkAny=True
                        regexPattern = ".*"+stemRegex+".*" 
                    elif curStem[0] == "-" == "-": # End
                        checkEnd=True
                        regexPattern = ".*"+stemRegex+"$"
                    elif curStem[len(curStem)-1] == "-": # Beg
                        checkBeg=True
                        regexPattern = "^"+stemRegex
                    
                    drug_match = drug.replace(" ","") 
                        
                    if (checkAny or checkBeg or checkEnd) and (re.search( regexPattern , drug_match)):
                        if len(result) == 0:
                            result.append({
                                "drug":drug,
                                "description":[row["definition"]],
                                "type":[curClassType]                                                           
                            })
                        else:
                            if result[len(result)-1]["drug"] == drug:
                                # verify that the description isn't here already.
                                definition_already_in = any(pd.DataFrame(result[len(result)-1])["description"].str.contains(row["definition"],regex=False))
                                if not definition_already_in:
                                    result[len(result)-1]["description"].append(row["definition"])
                                    result[len(result)-1]["type"].append(curClassType)
                            else:
                                result.append({
                                    "drug":drug,
                                    "description":[row["definition"]],
                                    "type":[curClassType]
                                })
        else:
            break
    return result

drug_classes = set_drugs_classes(l_task1 ,f_usan_stems)
exploded_drugs_classes = pd.DataFrame(drug_classes).apply(pd.Series.explode)


def aggregate_task2_enriched(l_task1 ,usans_stems):
    result = []
    cntrT = 0
    panda_t1 = pd.DataFrame(l_task1)
    alt_panda_t1 = panda_t1.explode("drugs")
    all_drugs = list(set(list(alt_panda_t1["drugs"])))
    all_drugs.sort()
    all_drugs = ' '.join(all_drugs).split()
        
    cntrLoop = 0
    for indexP1,rowP1 in alt_panda_t1.iterrows():
        drug = rowP1["drugs"]
        nct_id = rowP1["nct_id"]
        if cntrLoop < len(all_drugs):
            cntrLoop+=1
            for index, row in f_usan_stems.iterrows():
                curName = row["name"]
                curStem = row["stem"]
                isNanName = False
                isStr = type(row["name"]) == str
                if not isStr:
                    isNanName = np.isnan(row["name"])
                isSubGroup = False
                if not isNanName:
                    isSubGroup = "subgroup" in row["name"]
                isSubClass = isNanName or isSubGroup
                curClassType = "class"
                if isSubClass:
                    curClassType="subClass"

                regexPattern = ""
                if isinstance(curStem, str):
                    curStem = curStem.replace("=","")
                    stemRegex = curStem.replace("-","")
                    checkAny=False
                    checkBeg=False
                    checkEnd=False
                    if curStem[0] == "-" and curStem[len(curStem)-1] == "-": # Any
                        checkAny=True
                        regexPattern = ".*"+stemRegex+".*" # doubt about this... should I set .* at the end or just *?
                    elif curStem[0] == "-" == "-": # End
                        checkEnd=True
                        regexPattern = ".*"+stemRegex+"$"
                    elif curStem[len(curStem)-1] == "-": # Beg
                        checkBeg=True
                        regexPattern = "^"+stemRegex
                    
                    drug_match = drug.replace(" ","") 
                        
                    if (checkAny or checkBeg or checkEnd) and (re.search( regexPattern , drug_match)):
                        if len(result) == 0:
                            result.append({
                                "drug":drug,
                                "description":[row["definition"]],
                                "type":[curClassType],
                                "nct_id":[nct_id]
                            })
                        else:
                            if result[len(result)-1]["drug"] == drug:
                                definition_already_in = any(pd.DataFrame(result[len(result)-1])["description"].str.contains(row["definition"],regex=False))
                                if not definition_already_in:
                                    result[len(result)-1]["description"].append(row["definition"])
                                    result[len(result)-1]["type"].append(curClassType)
                                    result[len(result)-1]["nct_id"].append(nct_id)
                            else:
                                result.append({
                                    "drug":drug,
                                    "description":[row["definition"]],
                                    "type":[curClassType],
                                    "nct_id":[nct_id]
                                })
        else:
            break
    return result

ag_tsk2 = aggregate_task2_enriched(l_task1, f_usan_stems)
exploded_drugs_classes_enriched = pd.DataFrame(ag_tsk2).apply(pd.Series.explode)
ordered_exploded_drugs_classes_enriched = exploded_drugs_classes_enriched.sort_values(by="description")

# Probably doable without a loop, but we are in a rush
result_task3 = []
ordered_exploded_drugs_classes_enriched.sort_values(by="description")
for index,row in  ordered_exploded_drugs_classes_enriched.iterrows():
    curDrug = row["drug"]
    curType = row["type"]
    curDescription = row["description"]
    curNct_id = row["nct_id"]
    
    if(len(result_task3)==0):
        result_task3.append({
            "description":curDescription,
            "type":curType,
            "trials":[curNct_id]
        })
    else:
        if result_task3[len(result_task3)-1]["description"] == curDescription:
            result_task3[len(result_task3)-1]["trials"].append(curNct_id)
        else:
            result_task3.append({
                "description":curDescription,
                "type":curType,
                "trials":[curNct_id]
            })
    
# pd.DataFrame(result_task3)
# write the output result_task3
# ~~ write the result. 
pd_task_output_3 = pd.DataFrame(result_task3)

#create JSON file 
json_file3 = pd_task_output_3.to_json(orient='records') 

#export JSON file
with open('data/pd_task_output_3.json', 'w') as f:
    f.write(json_file3)

# ~~~~ for task 4
print("Task 4")
# ## We want to take back the previous output, and set a different structure for trials that used pairs of drugs. 
exploded_pd_task_output_3 = pd_task_output_3.apply(pd.Series.explode)
ordered_exploded_pd_task_output_3 = exploded_pd_task_output_3.sort_values(by="trials")
ordered_exploded_pd_task_output_3['Counts'] = ordered_exploded_pd_task_output_3.groupby(['description'])['type'].transform('count')
doublons_ordered_exploded_pd_task_output_3 = ordered_exploded_pd_task_output_3.loc[ordered_exploded_pd_task_output_3.duplicated(subset="trials",keep=False)]

# ## Misunderstood what was asked first. We aim to get, for each pair of classes (no subClasses) their names and sums of trials counts
# aggregate_class_pairs_trials_counts function expects a Pandas Data Frame
def aggregate_class_pairs_trials_counts( pd_task_output_3 ):
    output = []
    for index,row in pd_task_output_3.iterrows():
        # print("numbers of comparisons to make: ", len(pd_task_output_3)-index)
        nextInd = index+1
        while nextInd<len(pd_task_output_3):
            # verify if both objects are classes (we ignore subClasses)
            curClass = row["type"] == "class"
            nextClass = pd_task_output_3.iloc[nextInd,:]["type"] == "class"
            # if(index==5 and nextInd == 7):print("curClass: ",curClass,", nextClass: ",nextClass,', row["type"]: ',row["type"],', pd_task_output_3.iloc[nextInd,:]["type"]: ',pd_task_output_3.iloc[nextInd,:]["type"]) # testing
            # About trial count... Should we verify if they are the same? 
            # i.e. a trial using the two drugs be counted as two?
            if curClass and nextClass:
                output.append({
                    "description_1": row["description"],
                    "description_2": pd_task_output_3.iloc[nextInd,:]["description"],
                    "trial_count": len(row["trials"])+len(pd_task_output_3.iloc[nextInd,:]["trials"])
                })
            nextInd+=1
    # Sorting according to trial_count before returning
    ordered_output = pd.DataFrame(output).sort_values(by="trial_count")
    return ordered_output

pandas_output_4 = aggregate_class_pairs_trials_counts( pd_task_output_3 )        

# write the output result_task4
# ~~ write the result.
#create JSON file 
json_file4 = pandas_output_4.to_json(orient='records') 
#export JSON file
with open('data/pd_task_output_4.json', 'w') as f:
    f.write(json_file4)


print("Program completed")
