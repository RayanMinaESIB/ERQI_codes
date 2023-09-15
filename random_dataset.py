
import os
from numpy import round as npround
from numpy import abs as npabs
from numpy import log,sqrt,power,exp,random,sum,array,ones,concatenate,arange,set_printoptions,Inf,empty,min,max,average,argsort,std,mean
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

# =================================================    Compute h-index =====================================================
def comp_h_ind(citlist):
    citlist.sort(reverse=True)
    fh_index = 0
    cc = 0
    for elem in citlist:
        if elem>=cc+1:
            fh_index += 1
            cc += 1
        else:
            cc += 1
            break
    return fh_index
# ============================================================================================================================

# =================================================    Compute fh-index =====================================================
def comp_fh_ind(fcitlist):
    fcitlist2 = sorted(fcitlist,reverse=True)
    fh_index = 0
    cc = 0
    for elem in fcitlist2:
        if elem>=cc+1:
            fh_index += 1
            cc += 1
        else:
            cc += 1
            break
    return fh_index
# ============================================================================================================================

# ==================================================  Compute new distance ===================================================
def comp_dist(citlist,h_ind):
    tempoList = citlist.copy()
    inilen = len(tempoList)
    for i in range(inilen):
        if tempoList[i]<(h_ind+1):
            idx = i
            break
    disth = 0
    for i in range(idx,inilen):
        while tempoList[i]<(h_ind+1):
            tempoList[i] += 1
            disth += 1
            nh_ind = comp_h_ind(tempoList)
        if nh_ind==(h_ind+1):
            break
    counter = 1
    for i in range(inilen):
        if (tempoList[i]-(h_ind+1))==-1:
            counter += 1
    return disth/counter
# ============================================================================================================================

FolderName = 'C:\\' # Put a valid folder name on your machine

os.system("cls")
authnameList = list()
dfprof = pd.DataFrame()
authaff = 'random'
collist = ['Author','Affiliation','Papers','Co-authors','Citations','PAR','CAR','CPR','h','hf','d','he','ERQI','Citation-List','Authors-List']
set_printoptions(linewidth=Inf)

papmax = 200
citmax = 100
authmax = 10

nbcrang1 = 20
nbcrang2 = 40
nbcrang3 = citmax-nbcrang1-nbcrang2
prob1 = (1/40)*ones((nbcrang1))
prob2 = (1/114)*ones((nbcrang2))
prob3 = ((1-sum(prob1)-sum(prob2))/nbcrang3)*ones((nbcrang3))
probacit = concatenate([prob1,prob2,prob3])
# print('\n Citations Probability Distribution 20 60 100:',[round(100*sum(prob1),1),round(100*sum(prob2),1),round(100*sum(prob3),1)],'\n')

nbprang1 = 10
nbprang2 = 90
nbprang3 = papmax-nbprang1-nbprang2
prob1 = (1/80)*ones((nbprang1))
prob2 = (1/120)*ones((nbprang2))
prob3 = ((1-sum(prob1)-sum(prob2))/nbprang3)*ones((nbprang3))
probapap = concatenate([prob1,prob2,prob3])
# print('\n Papers Probability Distribution 10 100 200:',[round(100*sum(prob1),1),round(100*sum(prob2),1),round(100*sum(prob3),1)],'\n')

samplesize = 1000
for ii in range(samplesize):
    authnameList.append('Author-'+str(ii+1))
    citlistauth = list()
    nbpap = int(random.choice(array(arange(1,papmax+1)),size=1,p=probapap))
    citlist = sorted(random.choice(array(arange(citmax)),size=nbpap-1,p=probacit),reverse=True)
    citlist.append(0)
    authlist = random.randint(1,authmax,nbpap)
    nbcit = sum(citlist)
    nbauth = sum(authlist)
    par_ind = npround(nbpap/sqrt(nbauth),2)
    car_ind = npround(nbcit/nbauth,2)
    cpr_ind = npround(sqrt(nbcit/nbpap),2)
    h_ind = comp_h_ind(citlist)
    disth = npround(comp_dist(citlist,h_ind),2)

    if disth==100:
        he_ind = h_ind
    else:
        he_ind = npround(h_ind+exp(-disth),2)    

    fcitlist = npround(citlist/authlist,2)
    hf_ind = comp_fh_ind(fcitlist)
    erqi = npround(power(1+he_ind*hf_ind*car_ind*par_ind*cpr_ind,1/5),2)

    dfc = pd.DataFrame([[authnameList[ii],authaff,nbpap,nbauth,nbcit,par_ind,car_ind,cpr_ind,h_ind,hf_ind,disth,he_ind,erqi,citlist,authlist]],columns=collist,index=None)
    dfprof = pd.concat([dfprof,dfc],ignore_index=True)


corrmat = empty((5,5))
for jj in range(5):
    corrmat[jj,jj]=1
corrmat[0,1] = npabs(npround(pearsonr(dfprof['PAR'],dfprof['CAR'])[0],2))
corrmat[0,2] = npabs(npround(pearsonr(dfprof['PAR'],dfprof['CPR'])[0],2))
corrmat[0,3] = npabs(npround(pearsonr(dfprof['PAR'],dfprof['hf'])[0],2))
corrmat[0,4] = npabs(npround(pearsonr(dfprof['PAR'],dfprof['he'])[0],2))
corrmat[1,2] = npabs(npround(pearsonr(dfprof['CAR'],dfprof['CPR'])[0],2))
corrmat[1,3] = npabs(npround(pearsonr(dfprof['CAR'],dfprof['hf'])[0],2))
corrmat[1,4] = npabs(npround(pearsonr(dfprof['CAR'],dfprof['he'])[0],2))
corrmat[2,3] = npabs(npround(pearsonr(dfprof['CPR'],dfprof['hf'])[0],2))
corrmat[2,4] = npabs(npround(pearsonr(dfprof['CPR'],dfprof['he'])[0],2))
corrmat[3,4] = npabs(npround(pearsonr(dfprof['hf'],dfprof['he'])[0],2))

for ii in range(5):
    for jj in range(ii,5):
        corrmat[jj,ii] = corrmat[ii,jj]
print('\n\n Correlation Matrix \n\n',corrmat)

alldata = dfprof.to_numpy()
alldata  = alldata[:,5:13]
print('\n============================================\n')

sorted_hidx = argsort(-alldata[:,3])
alldata_sorted_h = alldata[sorted_hidx]

sorted_erqi = argsort(-alldata[:,7])
alldata_sorted_erqi = alldata[sorted_erqi]
h_data = alldata_sorted_h[:,3]
erqi_data = alldata_sorted_erqi[:,7]

hrank = arange(1,samplesize+1).reshape((samplesize,1))
sorted_h = concatenate((hrank,alldata_sorted_h),axis=1)

erqiidx = argsort(-sorted_h[:,8])
sorted_erqi_rank = sorted_h[erqiidx]

hrank_univ1 = arange(1,31+1).reshape((31,1))
sorted_erqi_rank_univ1 = [1,2,3,8,5,22,15,4,6,10,9,18,7,12,14,11,16,17,19,23,13,20,21,27,24,25,26,28,30,29,31]
erqi_data_univ1 = [9.32,5.94,4.91,4.53,4.47,3.98,3.96,3.84,3.77,3.46,3.39,3.38,3.36,3.33,3.23,3.17,3.14,2.91,2.84,2.49,2.42,2.33,2.25,1.60,1.59,1.56,1.48,1.20,1.15,1.00,1.00]
h_data_univ1 = [29,17,15,11,11,10,9,9,8,8,7,6,6,6,6,6,6,6,5,5,5,4,4,4,3,2,2,2,1,1,0]

fsleg = 14
fstit = 13
fslab = 13
plt.figure(figsize=(8,6),dpi=400)  # figsize=(10,6),dpi=400
plt.subplot(1,2,1)
plt.plot(hrank_univ1,sorted_erqi_rank_univ1,'-ro',linewidth=2.5)
plt.plot(hrank_univ1,hrank_univ1,'-b',linewidth=3)
plt.xlabel('ERQI ranks',fontsize=fslab)
plt.grid(visible=True,which='major',axis='y',color='#1F497D',linestyle=':')
plt.grid(visible=True,which='major',axis='x',color='#1F497D',linestyle=':')
plt.axis([0,31,1,31])
plt.legend(['h ranking','ERQI ranking'],loc='best',fontsize=fsleg)
plt.title('Univ-1 Dataset',fontsize=fstit,color='#1F497D')
plt.xticks(fontsize=fslab)
plt.yticks(fontsize=fslab)

plt.subplot(1,2,2)
plt.plot(hrank,sorted_erqi_rank[:,0],'-r',linewidth=2.5)
plt.plot(hrank,hrank,'-b',linewidth=3)
plt.xlabel('ERQI ranks',fontsize=fslab)
plt.grid(visible=True,which='major',axis='y',color='#1F497D',linestyle=':')
plt.grid(visible=True,which='major',axis='x',color='#1F497D',linestyle=':')
plt.axis([0,samplesize,1,samplesize])
plt.legend(['h ranking','ERQI ranking'],loc='best',fontsize=fsleg)
plt.title('Randomly Generated Dataset',fontsize=fstit,color='#1F497D')
plt.xticks(fontsize=fslab)
plt.yticks(fontsize=fslab)
plt.savefig(os.path.join(FolderName,'Figure_4.png'))

plt.figure(figsize=(8,6),dpi=400)
plt.subplot(1,2,1)
plt.plot(hrank_univ1,h_data_univ1,'-ro',linewidth=2,markerfacecolor='None')
plt.plot(hrank_univ1,erqi_data_univ1,'-bo',linewidth=2,markerfacecolor='None')
plt.xlabel('ERQI ranks',fontsize=fslab)
plt.grid(visible=True,which='major',axis='y',color='#1F497D',linestyle=':')
plt.grid(visible=True,which='major',axis='x',color='#1F497D',linestyle=':')
plt.axis([0,31,1,30])
plt.legend(['h values','ERQI values'],loc='best',fontsize=fsleg)
plt.title('Univ-1 Dataset',fontsize=fstit,color='#1F497D')
plt.xticks(fontsize=fslab)
plt.yticks(fontsize=fslab)

plt.subplot(1,2,2)
plt.plot(hrank,h_data,'-r',linewidth=2.5)
plt.plot(hrank,erqi_data,'-b',linewidth=3)
plt.xlabel('ERQI ranks',fontsize=fslab)
plt.grid(visible=True,which='major',axis='y',color='#1F497D',linestyle=':')
plt.grid(visible=True,which='major',axis='x',color='#1F497D',linestyle=':')
plt.axis([0,samplesize,1,50])
plt.legend(['h values','ERQI values'],loc='best',fontsize=fsleg)
plt.title('Randomly Generated Dataset',fontsize=fstit,color='#1F497D')
plt.xticks(fontsize=fslab)
plt.yticks(fontsize=fslab)
plt.savefig(os.path.join(FolderName,'Figure_2.png'))

print('\n h index: mu and sigma --> ',npround(mean(h_data),1),npround(std(h_data),1))
print('\n ERQI: mu and sigma --> ',npround(mean(erqi_data),1),npround(std(erqi_data),1))


plt.figure(figsize=(5,4),dpi=300)
plt.boxplot(h_data,positions=[1],patch_artist=True,boxprops=dict(facecolor='red',color='black'),whiskerprops=dict(color='red'),capprops=dict(color='red'),medianprops=dict(color='yellow'))
plt.grid(visible=True,which='major',axis='y',color='#1F497D',linestyle=':')
plt.grid(visible=True,which='major',axis='x',color='#1F497D',linestyle=':')
plt.boxplot(erqi_data,positions=[2],patch_artist=True,boxprops=dict(facecolor='blue',color='black'),whiskerprops=dict(color='blue'),capprops=dict(color='blue'),medianprops=dict(color='yellow'))
plt.title('Boxplots: h index (1,red) and ERQI (2,blue)')

plt.savefig(os.path.join(FolderName,'Figure_3.png'))


writer = pd.ExcelWriter(os.path.join(FolderName,'Random_Dataset.xlsx'))
dfprof.to_excel(writer,sheet_name='Sheet1')
workbook  = writer.book
worksheet = writer.sheets['Sheet1']
worksheet.freeze_panes(1,0)
worksheet.set_column('A:A',5)
worksheet.set_column('B:B',20)
worksheet.set_column('C:F',10)
worksheet.set_column('G:N',8)
worksheet.set_column('O:O',100)
data_format1 = workbook.add_format({'bg_color': '#DDF0C8','reading_order':'2','align':'center','border': 1,'font_size':13})
data_format2 = workbook.add_format({'bg_color': '#FFFFFF','reading_order':'2','align':'center','border': 1,'font_size':13})
for row in range(1,dfprof.shape[0]+1,2):
    worksheet.set_row(row,cell_format=data_format1)
for row in range(2,dfprof.shape[0]+1,2):
    worksheet.set_row(row,cell_format=data_format2)
writer.save()
