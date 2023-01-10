import pandas as pd
from ms2pip.single_prediction import SinglePrediction
import spectrum_utils.spectrum as sus
import math
from itertools import permutations
from copy import deepcopy

def switch_pos_center(seq,pos_switch="middle",rev_length=2):
    if pos_switch == "middle":
        sel_pos = int(len(seq)/2)
    if pos_switch == "left":
        sel_pos = int(int(len(seq)/2)/2)
    if pos_switch == "right":
        sel_pos = int(len(seq)/2)+int(int(len(seq)/2)/2)
    new_seq = seq[:sel_pos-math.ceil(rev_length/2)+1]+seq[sel_pos-math.ceil(rev_length/2)+1:sel_pos+int(rev_length/2)+1][::-1]+seq[sel_pos+int(rev_length/2)+1:]
    return pd.Series({"seq":new_seq,"aa_sw1":seq[sel_pos-math.ceil(rev_length/2)+1:sel_pos+int(rev_length/2)+1],"aa_sw2":seq[sel_pos-math.ceil(rev_length/2)+1:sel_pos+int(rev_length/2)+1][::-1]})

def get_all_combs(seq,rev_length=2,pos_switch="middle"):
    if pos_switch == "middle":
        sel_pos = int(len(seq)/2)
    if pos_switch == "left":
        sel_pos = int(int(len(seq)/2)/2)
    if pos_switch == "right":
        sel_pos = int(len(seq)/2)+int(int(len(seq)/2)/2)
    
    sel_seq = seq[sel_pos-math.ceil(rev_length/2)+1:sel_pos+int(rev_length/2)+1][::-1]
    all_combs = ["".join(v) for v in permutations(sel_seq)]

    ret_dict = {}
    for idx,v in enumerate(all_combs):
        new_seq = seq[:sel_pos-math.ceil(rev_length/2)+1]+v+seq[sel_pos+int(rev_length/2)+1:]
        if v == sel_seq:
            continue
        ret_dict[idx] = {"seq":new_seq,"aa_sw1":sel_seq,"aa_sw2":v}

    return pd.DataFrame(ret_dict)


def mutate_aa(seq,mutate_to):
    sel_pos = int(len(seq)/2)
    orig_aa = str(seq[sel_pos])
    new_seq = list(seq)

    new_seq[sel_pos] = mutate_to
    new_seq = "".join(new_seq)
    return pd.Series({"seq":new_seq,"sel_pos":sel_pos,"aa_sw1":orig_aa,"aa_sw2":new_seq[sel_pos]})

#def switch_two_pos_center(seq,pos_switch="middle"):
#    sel_pos = int(len(seq)/2)
#    return pd.Series({"seq":seq[:sel_pos]+seq[sel_pos:sel_pos+4][::-1]+seq[sel_pos+4:],"aa_sw1":}

def get_predicted_spectrum(peptide, modifications, charge, model="HCD2021", ms2pip_instance=None):
    if not ms2pip_instance:
        ms2pip = SinglePrediction()
    else:
        ms2pip = ms2pip_instance
        
    mz, intensity, annotation = ms2pip.predict(peptide, modifications, charge, model=model)
    
    return (mz, intensity)