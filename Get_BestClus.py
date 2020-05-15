#!/usr/bin/env python
import sys

def res_dict(dict_file):
	'''This function uses a .txt file that contains the PDBID's in the first column, and the resolution in the second column and converts it into a dictionary''' 
	with open(dict_file) as res:
		res_dic = {}
		header = res.readline()
		for i in res:
			i = i.rstrip().rsplit()
			res_dic[i[0]]=float(i[1])
	return(res_dic)


def min_ids(clus_file,res_dic):
	'''this function takes a file containing the cluster members line by line, separated by a tab and returns the element of each cluster with the highest resolution (lowest value)'''
	with open (clus_file) as cluster:
		min_id = []
		for i in cluster:
			i = i.rstrip().rsplit('\t') 
			#print(i)
			min = 3.5
			min_id = ''
			for j in range(len(i)):
				if res_dic[i[j]] < min:
					min = res_dic[i[j]]
					min_id = i[j] 
			print(min_id,":",min)

if __name__ == '__main__':
	filename1 = sys.argv[1]
	filename2 = sys.argv[2]
	print(min_ids(filename1,res_dict(filename2)))