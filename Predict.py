# -*- coding: cp936 -*-
import xlrd
#import xlwt
from numpy import *
import datetime
import os

#initialize parameter:
starttime=datetime.datetime.now()
Similarity_threshold=0.2    #smaller than 0.5,don't show the links
Pa=2.26
Max_length=3       #Maximun path length
index1=0

#input the name lists of diseases and lncRNAs
print "Step 1. input the name lists of diseases and lncRNAs.."
data = xlrd.open_workbook("disease_name.xlsx")
table = data.sheets()[0]
x = table.nrows
disease_data = []
for index in range(x):
    disease_data.append([index+1,table.row(index)[0].value.encode("UTF-8")])

data = xlrd.open_workbook("lncRNA_name.xlsx")
table = data.sheets()[0]
y = table.nrows
lncRNA_data = []
for index in range(y):
    lncRNA_data.append([index+1,table.row(index)[0].value.encode("UTF-8")])


print "Step 2. input other information resource.."
DSM=zeros([2,x,x])  #disease similarity matrices
LSM=zeros([2,y,y])  #lncRNA similarity matrices
DSM1=zeros([x,x]) #1 OR 0
LSM1=zeros([y,y])

#judge whether d(i) and d(j) has disease semantic similarity
data = xlrd.open_workbook("Judgedd.xlsx")
table = data.sheets()[0]
for i in range(table.nrows):
    for j in range(table.ncols):
        DSM1[i][j] = table.row(i)[j].value

#judge whether l(i) and l(j) has miRNA functional similarity
data = xlrd.open_workbook("Judgell.xlsx")
table = data.sheets()[0]
for i in range(table.nrows):
    for j in range(table.ncols):
        LSM1[i][j] = table.row(i)[j].value

#Type 1 disease semantic similarity
data = xlrd.open_workbook("DD.xlsx")
table = data.sheets()[0]
for i in range(table.nrows):
    for j in range(table.ncols):
        DSM[0][i][j] = table.row(i)[j].value

#lncRNA functional similarity
data = xlrd.open_workbook("LL.xlsx")
table = data.sheets()[0]
for i in range(table.nrows):
    for j in range(table.ncols):
        LSM[0][i][j] = table.row(i)[j].value


while index1 < 1:
    print "Step 3. complete the construction of network by the integration"
    
    data = xlrd.open_workbook("Gaussion_disease.xlsx")
    table = data.sheets()[0]
    for i in range(table.nrows):
        for j in range(table.ncols):
            DSM[1][i][j] = table.row(i)[j].value
    
    data = xlrd.open_workbook("Gaussion_lncRNA.xlsx")
    table = data.sheets()[0]
    for i in range(table.nrows):
        for j in range(table.ncols):
            LSM[1][i][j] = table.row(i)[j].value
    
    A = zeros([y,y])
    B = zeros([x,x])
    
    for i in range(y):
        for j in range(y):
            if LSM1[i][j] == 1:
                A[i][j] = LSM[0][i][j]
            else:
                A[i][j] = LSM[1][i][j]
    
    for i in range(x):
        for j in range(x):
            if DSM1[i][j] == 1:
                B[i][j] = DSM[0][i][j]
            else:
                B[i][j] = DSM[1][i][j]
    
    "print Step 4. build the adajency matrix Y.."
    data = xlrd.open_workbook("LD.xlsx")
    table = data.sheets()[0]
    Y = zeros([y,x])
    for i in range(table.nrows):
        for j in range(table.ncols):
            Y[i][j] = table.row(i)[j].value
    
    
    print "Step 5. constructed a heterogeneous graph consisting of three interlinked sub-graphs.."
    Network_disease = zeros([x,4,max(x,y)])
    Network_lncRNA = zeros([y,4,max(x,y)])
    
    for i in range(x):
        index = 0
        for j in range(x):
            if B[i][j] > Similarity_threshold:
                Network_disease[i][0][index] = B[i][j]
                Network_disease[i][1][index] = j+1
                index += 1
        index = 0 
        for k in range(y):
            if Y[k][i] == 1:
                Network_disease[i][2][index] = 1
                Network_disease[i][3][index] = k+1
                index += 1
    for i in range(y):
        index = 0
        for j in range(y):
            if A[i][j] > Similarity_threshold:
                Network_lncRNA[i][2][index] = A[i][j]
                Network_lncRNA[i][3][index] = j+1
                index += 1
        index = 0
        for k in range(x):
            if Y[i][k] == 1:
                Network_lncRNA[i][0][index] = 1
                Network_lncRNA[i][1][index] = k+1
                index += 1
    
    print "Step 6. start the depth-first search algorithm.."
    P = zeros([x,y])
    for i in range(x):
        ll = []
        lll = [i+1]
        
        for k in range(max(y,x)): #disease
            if Network_disease[i][0][k] == 0 :
                break
            if int(Network_disease[i][1][k]) != i+1:
                lll.append(int(Network_disease[i][1][k]))
                ll.append(lll)
                lll = [i+1]
        for k in range(max(y,x)): #lncRNA
            if Network_disease[i][2][k] == 0:
                break
            lll.append(-int(Network_disease[i][3][k]))
            ll.append(lll)
            lll = [i+1]
            P[i][-ll[-1][-1]-1] += 1
        
        for j in range(Max_length-1):
            ll1 = []
            for k in range(len(ll)):
                
                lll = ll[k]
                if ll[k][j+1] > 0 : #disease
                    for m in range(max(y,x)): #disease
                            if Network_disease[ll[k][j+1]-1][0][m] == 0:
                                break;
                            if (int(Network_disease[ll[k][j+1]-1][1][m])) not in ll[k] :
                                ll1.append(ll[k]+[int(Network_disease[ll[k][j+1]-1][1][m])])
                    for m in range(max(y,x)): #lncRNA
                            temp = 1
                            if Network_disease[ll[k][j+1]-1][2][m] == 0:
                                break;
                            if (-int(Network_disease[ll[k][j+1]-1][3][m])) not in ll[k] :
                                ll1.append(ll[k]+[-int(Network_disease[ll[k][j+1]-1][3][m])])
                                for q in range(j+2):
                                    if ll1[-1][q+1] > 0 and ll1[-1][q] > 0: #disease-disease
                                        temp *= B[ll1[-1][q]-1][ll1[-1][q+1]-1]
                                    elif ll1[-1][q+1] < 0 and ll1[-1][q] > 0: #disease-lncRNA
                                        temp *= 1
                                    elif ll1[-1][q+1] > 0 and ll1[-1][q] < 0: #lncRNA-disease
                                        temp *= 1
                                    elif ll1[-1][q+1] < 0 and ll1[-1][q] < 0: #lncRNA-lncRNA
                                        temp *= A[-ll1[-1][q]-1][-ll1[-1][q+1]-1]
                                P[i][-ll1[-1][-1]-1] += (temp**(Pa*(j+2)))
            
    
                if ll[k][j+1] < 0: #lncRNA
                    for m in range(max(y,x)): #disease
                            if Network_lncRNA[-ll[k][j+1]-1][0][m] == 0:
                                break;
                            if (int(Network_lncRNA[-ll[k][j+1]-1][1][m])) not in ll[k] :
                                ll1.append(ll[k]+[int(Network_lncRNA[-ll[k][j+1]-1][1][m])])
                    for m in range(max(y,x)): #lncRNA
                            temp = 1
                            if Network_lncRNA[-ll[k][j+1]-1][2][m] == 0:
                                break;
                            if (-int(Network_lncRNA[-ll[k][j+1]-1][3][m])) not in ll[k]:
                                ll1.append(ll[k]+[-int(Network_lncRNA[-ll[k][j+1]-1][3][m])])
                                for q in range(j+2):
                                    if ll1[-1][q+1] > 0 and ll1[-1][q] > 0: #disease-disease
                                        temp *= B[ll1[-1][q]-1][ll1[-1][q+1]-1]
                                    elif ll1[-1][q+1] < 0 and ll1[-1][q] > 0: #disease-lncRNA
                                        temp *= 1
                                    elif ll1[-1][q+1] > 0 and ll1[-1][q] < 0: #lncRNA-disease
                                        temp *= 1
                                    elif ll1[-1][q+1] < 0 and ll1[-1][q] < 0: #lncRNA-lncRNA
                                        temp *= A[-ll1[-1][q]-1][-ll1[-1][q+1]-1]
                                P[i][-ll1[-1][-1]-1] += (temp**(Pa*(j+2)))
            
            ll = ll1         
    #f = xlwt.Workbook()
    #sheet = f.add_sheet("sheet1")
    #for i in range(y):
        #for j in range(x):
            #sheet.write(i,j,P[j][i])
    #f.save("predictLDABPL_result2.xls")
    
    print "Step 7. construct the dictionary for ranking.."
    dic=[]
    for i in range(y):
        for j in range(x):
            if Y[i][j] != 1:
            #get rid of the known lncRNA-disease associations, only prioritize the unverified lncRNA-disease pairs
            #P matrix records the aggregated score for each lncRNA-disease pair
                unit=tuple((lncRNA_data[i][1],disease_data[j][1],P[j][i],lncRNA_data[i][0],disease_data[j][0]))
                dic.append(unit)
    #Rank the unverified lncRNA-disease pairs based on their predicted scores            
    dic=sorted(dic,key=lambda x:x[2])
    dic.reverse()
    
####################################################################    
 

####################################################################
#Store the global ranking result for all given diseases (disease_name, lncRNA_name, score)
    print "Step 8. store the global ranking result.."
    f0=open("Ranking_result.txt",'w')
    
    for i in range(len(dic)):
        f0.write(str(dic[i][1])+"\t"+str(dic[i][0])+"\t"+str(dic[i][2])+"\n")
        f0.write('\n')
    f0.close()
####################################################################

####################################################################
#Store the ranking result for each diseases in the local dictory "Prediction"
#name each profile with the disease name (lncRNA_name, score)
    print "Step 9. store the each disease profile.."
    os.mkdir("Prediction")
    for i in range(len(dic)):
        f=open("Prediction/"+str(dic[i][1])+".txt",'a')
        f.write(str(dic[i][1])+"\t"+str(dic[i][0])+"\t"+str(dic[i][2])+"\n")
        f.close()
        
    index1+=1
    endtime=datetime.datetime.now()
    print "Total runtime "+str((endtime-starttime).seconds)+" second"
print "Completed!"        
####################################################################