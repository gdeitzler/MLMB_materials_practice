
"""
@author: Tatiana Lenskaia

Updated: May 7, 2021

Python3
"""

import random;



def MFastaCutter(fInName):
    #Multi fasta cutter
    
    t_names = []
    text = ""
    name = ""
    
    fIn = open(fInName,"r")
    for line in fIn:
        if ">" in line:
            if text != "" and name != "":
                fOut = open(name+".fasta","w")
                fOut.write(text)
                fOut.close()
            text = line
            name = line[1:].split(" ")[0]
            t_names.append(name)
        else:
            text = text+line

    
    if text != "" and name != "":
        fOut = open(name+".fasta","w")
        fOut.write(text)
        fOut.close()
    
    fname = fInName.rsplit(".",1)[0]+".txt"
    fOut = open(fname, "w")
    s = ""
    for it in t_names:
        s = s+it+"\n"
    fOut.write(s[:(-1)])
    fOut.close()
    return fname


def CountGC(seq):
    seq = seq.lower();
    n_seq = len(seq)
    
    n_a = seq.count("a");
    n_c = seq.count("c");
    n_g = seq.count("g");
    n_t = seq.count("t");
    
    n_bases = n_a + n_c + n_g + n_t
    
    n_other = n_seq - n_bases;
    
    if n_bases != 0:
        gc = round(1.0*(n_c+n_g)/n_bases*100,2)
    else:
        gc = "N/A"
    
    return [gc, n_other]



def GetListFromFile(fInName):
    t = []
    fIn = open(fInName, "r")
    for line in fIn:
        line = line.strip()
        if line != "":
            if line not in t:
                t.append(line)
    fIn.close()
    return t



#Updated: May 4, 2019   
def GetText(finName):
    '''Extracts text from a single fasta file'''
    fin = open(finName, 'r')
    text = ''
    for line in fin:
        line = line.strip()
        if (line != "") and (line[0] != '>'):
            text = text + line
    fin.close()
    return text


def CreateDict_count(text, m, gtp = "l"):
    if (m <= 0) or (m > len(text)):
        print("(n = "+str(m)+") is not a valid window size for this genome!!!");
        return {}
	
    d_g = dict()
    nn = len(text);
    gtype = gtp.lower();
    gtype = gtype[0];
	
    if gtype == "c":
        text = text + text[0:(m-1)];
        lastpos = nn;
    elif gtype == "l":
        lastpos = nn-m+1;
    else:
        print("Is this genome linear or circular?");
        return d_g;
		
    for ii in range (lastpos):
        bl = text[ii:(ii+m)]; 
        bl = bl.upper();
        if bl in d_g :
            d_g[bl] = d_g[bl] + 1
        else:
            d_g[bl] = 1
    return d_g;	




def FindIntersection(d_g11, d_g22):
	
	t_ints = list()
	nn1 = len(d_g11);
	nn2 = len(d_g22);
	
	if nn1 <= nn2:
		d_first = d_g11;
		d_last = d_g22;
	else:
		d_first = d_g22;
		d_last = d_g11;
	
	for s in d_first:
		if s in d_last:
			t_ints.append(s);
	return t_ints;


    
    
    
 
    
def PrintMatrix(k, pref, row_index, col_index, matrix, sep = ","):
    fOut = open(str(k)+pref+".csv","w");
    fOut.write(str(k)+sep+sep.join(str(x) for x in col_index)+"\n")
    for i in range(len(matrix)):
        s = str(row_index[i])
        for el in matrix[i]:
            s = s+sep+str(el)
        fOut.write(s+"\n")
        
    fOut.close()
    return

