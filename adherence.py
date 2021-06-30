# ADHERENCE.PY - API for generating adherence metrics
# Author: Thomas Patton
#
# (c) 2021 Cohere Health
#
#
# How to use:
# 1. Call the initialize() function. This will pull the appropriate files from S3.
# 2. Determine which of the adherence functions is right for your use case
# 3. Use the create_udf() function on the chosen adherence function. This will return a PySpark UDF which
#    can be applied to your claims dataframe.

import os
import boto3
import matplotlib.pyplot as plt
import pandas
import numpy
from io import StringIO
from Bio import pairwise2
from pyspark.sql.types import FloatType, DoubleType, StringType
from pyspark.sql import functions as F
import ast
import json


def initialize():
    '''
    Load resources from S3
    '''
    s3 = boto3.resource('s3')
    
    json_str = s3.Object('sagemaker-coherehealth-1', 'claims-data/carepaths.json').get()['Body'].read().decode('utf-8')
    data = json.loads(json_str)

    hierarchy_str = s3.Object('sagemaker-coherehealth-1', 'claims-data/hierarchy.csv').get()['Body'].read().decode('utf-8')
    hdf = pd.read_csv(StringIO(hierarchy_str))
    hdf.at[112, 'End Code'] = 99499
    hdf = hdf.astype({'End Code': 'int64'})
    hierarchy = list(hdf.T.to_dict().values())

    cat_map = get_category_map(hierarchy)
    

def get_category_map(hierarchy):
    '''
    Creates a dictionary to map CPT codes to minor categories
    '''
    cat_map = {}
    for h in hierarchy:
        for i in range(h['Start Code'], h['End Code'] + 1):
            cat_map[str(i)] = h['CPT Minor Category']
    return cat_map


def make_scoring_dict(target_str):
    '''
    Creates a scoring dictionary for the sequence alignment algorithm
    '''
    sdict = {}
    for c1 in target_str:
        for c2 in target_str:
            if c1 == c2:
                if c1 != '0' and c2 != '0':
                    sdict[(c1, c2)] = 10
                else:
                    sdict[(c1, c2)] = 1000
            elif c1 == '0' or c2 == '0':
                sdict[(c1, c2)] = -1000
            elif abs(int(c1) - int(c2)) == 1:
                sdict[(c1, c2)] = 5
            else:
                sdict[(c1, c2)] = 0
    return sdict


def get_target_sequence(carepath):
    '''
    Generates the "target" sequence
    '''
    mapping = {}
    key_list = list(carepath.keys())
    target = '0'
    for k in range(len(key_list)):
        mapping[key_list[k]] = str(k+1)
        target += str(k+1)
        
    return target, mapping


def get_reduced_str(input_str):
    '''
    Reduce string to only unique characters
    '''
    new = ''
    for c in input_str:
        if c not in new:
            new += c
    return new


def construct_pcd_seq(procedures, carepath, mapping):
    '''
    function to query Px codes and get their service cateogry / construct Px sequence
    '''
    pcd_seq = '0'
    for procedure in procedures:
        for key in carepath:
            if procedure in carepath[key]:
                pcd_seq += mapping[key]
                break
                
    return pcd_seq, (len(pcd_seq)-1) / len(procedures)


def code_to_category(code, hierarchy):
    '''
    Gets the minor category associated with a code
    '''
    for h in hierarchy:
        if code >= h['Start Code'] and code <= h['End Code']:
            return h['CPT Minor Category']
    return 'unknown'


def get_category_map(hierarchy):
    '''
    Creates a dictionary mapping CPT codes to their minor categories
    '''
    cat_map = {}
    for h in hierarchy:
        for i in range(h['Start Code'], h['End Code'] + 1):
            cat_map[str(i)] = h['CPT Minor Category']
    return cat_map


def construct_rollup_pcd_seq(procedures, carepath, mapping, cat_map):
    '''
    Attmept to find an exact match. If none exists, attempt to roll up
    '''
    pcd_seq = '0'
    for procedure in procedures:
        if procedure == 'None' or procedure == '':
            continue
            
        exact_match = False
        for key in carepath:
            if procedure in carepath[key]:
                pcd_seq += mapping[key]
                exact_match = True
                break
        
        # no direct match, attempt to roll up
        if not exact_match and procedure.isdigit():
            rollup_match = False
            try:
                cat = cat_map[procedure.lstrip('0')]
            except:
                continue
            for key in carepath:
                for cp_prcd in carepath[key]:
                    if cp_prcd.isdigit():
                        try:
                            prcd_cat = cat_map[cp_prcd.lstrip('0')]
                        except:
                            continue
                        if cat is cat_map[cp_prcd.lstrip('0')]:
                            pcd_seq += mapping[key]
                            rollup_match = True
                            break
                if rollup_match:
                    break
    return pcd_seq, (len(pcd_seq)-1) / len(procedures)


def sequence_alignment(pcd_seq, cp_seq):
    '''
    Modified sequence alignment
    '''
    scoring = make_scoring_dict(cp_seq)
    a = pairwise2.align.globaldd(cp_seq, 
                                 pcd_seq,
                                 scoring,
                                 -10,
                                 0,
                                 -5,
                                 -1,
                                 penalize_end_gaps=False,
                                 score_only=True)
    
    return ((a-1000) / (len(pcd_seq)-1) / 10)


def seq_alignment_adherence(diag_cd, procedure_cds):
    try:
        carepaths = data[diag_cd][0]
        best_adh = 0.0
        for cp in carepaths:
            target, mapping = get_target_sequence(cp)
            code_list = procedure_cds.split(' ')[:-1]
            pcd_seq, frac = construct_pcd_seq(code_list, cp, mapping)
            reduced_pcd_seq = get_reduced_str(pcd_seq)
            if frac != 0:
                adherence = sequence_alignment(reduced_pcd_seq, target) * frac
                print(adherence)
            else:
                adherence = 0.0
            if adherence > best_adh:
                best_adh = adherence
        return round(best_adh, 3)
    except:
        return -1.0

    
def generous_adherence(diag_cd, procedure_cds):
    try:
        carepaths = data[diag_cd][0]
        best_adh = 0.0
        for cp in carepaths:
            target, mapping = get_target_sequence(cp)
            code_list = procedure_cds.split(' ')[:-1]
            pcd_seq, frac = construct_pcd_seq(code_list, cp, mapping)
            reduced_pcd_seq = get_reduced_str(pcd_seq)
            if frac != 0:
                adherence = sequence_alignment(reduced_pcd_seq, target)
            else:
                best_adh = -1.0
            if adherence > best_adh:
                best_adh = adherence
        return round(best_adh, 3)
    except:
        return -1.0
        
    
def rollup_adherence(diag_cd, procedure_cds):
    try:
        carepaths = data[diag_cd][0]
        best_adh = 0.0
        for cp in carepaths:
            target, mapping = get_target_sequence(cp)
            code_list = procedure_cds.split(' ')[:-1]
            pcd_seq, frac = construct_rollup_pcd_seq2(code_list, cp, mapping, cat_map)
            reduced_pcd_seq = get_reduced_str(pcd_seq)
            if frac != 0:
                adherence = sequence_alignment(reduced_pcd_seq, target) * frac
            else:
                adherence = 0.0
            if adherence > best_adh:
                best_adh = adherence
        return round(best_adh, 3)
    except Exception as e:
        return -1.0
    
    
def create_udf(fn):
    '''
    Creates a UDF of a given rollup function
    '''
    return F.udf(fn, DoubleType())