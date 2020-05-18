#!/usr/bin/python
import sys
import math

def get_hmmout(filename):
    '''Takes as input a file containing 3 columns: col 1 contains the Uniprot ID's, col 2 contains the E-value (either sequence or domain E-values resulting from a hmmsearch), 
    and the col 3 a '0' for the negative set, and a '1' for the negative set'''
    f_list=[]
    d={}
    with open (filename) as f1:
        for i in f1:
            v=i.rstrip().split()
            d[v[0]]=d.get(v[0],[])
            d[v[0]].append([float(v[1]),int(v[2])])
        for ids in d.keys():
            #we may get several hits for one single UniProtKb Id, then we only want to keep the one with the lowest E-value.
            d[ids].sort()
            f_list.append(d[ids][0])
        return(f_list)

def get_data(filename):
    ldata=[]
    with open(filename) as f2:
        for i in f2:
            v=i.rstrip().split()
            ldata.append([float(v[1]),int(v[2])])
        return(ldata)

def get_coma(data,t):
    '''from a file containing the e-values and the model we can calculate the confusion matrix
    t = the selected threshold to calculate the values that are above and below the threshold'''
    cm = [[0.0,0.0],[0.0,0.0]]
    #we build a confusion matrix with starting values 0 
    for i in data:
        if i[0]<t and i[1]==1: #if e-value < treshold and belongs to the positive set
            cm[0][0]=cm[0][0] + 1 #cm[0][0] +1 (true positives: TP)
        if i[0]>=t and i[1]==1: #if e-value ≥ treshold and positive set
            cm[1][0]=cm[1][0] + 1 #cm[1][0] +1 (false negatives: FN)
        if i[0]<t and i[1]==0: #if e-value < treshold and belongs to the negative set
            cm[0][1]=cm[0][1]+1 #cm[0][1] +1 (false positives: FP)
        if i[0]>=t and i[1] ==0: #if e-value > treshold and belongs to the negative set
            cm[1][1] = cm[1][1]+1 #cm[1][1] + 1 (true negatives: TN)
    return(cm)

def get_ac(m):
    '''it takes a confusion matrix and calculates the accuracy'''
    return(float(m[0][0]+m[1][1])/(sum(m[0])+sum(m[1])))
    #(TP + TN) / (TP + FP + FN + TN) 

def mcc(n):
    '''takes a confusion matrix and calculates the Matthew Correlation Coefficient'''
  d=(n[0][0]+n[1][0])*(n[0][0]+n[0][1])*(n[1][1]+n[1][0])*(n[1][1]+n[0][1])
  return (n[0][0]*n[1][1]-n[0][1]*n[1][0])/math.sqrt(d)

if __name__== "__main__":
    filename=sys.argv[1]
    #th=float(sys.argv[2])
    #sel_t =  10**-9 #this field can be modified with a specific treshold to calculate acc and mcc from one of the sets
    data = get_hmmout(filename)
    for i in range(20): #iterates to 20 different e-values
        t=10**-i
        cm= get_coma(data,t)
        print ('TH:',t,'ACC:',get_ac(cm),'MCC:',mcc(cm), cm)
    #cm = get_coma(data,10**-9)
    #print('TH:',sel_t,'ACC:',get_ac(cm),'MCC:',mcc(cm), cm)