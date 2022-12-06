import pandas as pd

def switch_two_pos(seq,pos_switch="middle"):
    sel_pos = int(len(seq)/2)
    return seq[:sel_pos]+seq[sel_pos:sel_pos+2][::-1]+seq[sel_pos+2:]

def switch_two_pos_center(seq,pos_switch="middle"):
    sel_pos = int(len(seq)/2)
    return seq[:sel_pos]+seq[sel_pos:sel_pos+3][::-1]+seq[sel_pos+3:]

def switch_two_pos_center(seq,pos_switch="middle"):
    sel_pos = int(len(seq)/2)
    return seq[:sel_pos]+seq[sel_pos:sel_pos+4][::-1]+seq[sel_pos+4:]