import math


list_GBDT = [ 1.121, 1.370, 1.394, 1.494]
list_QBDT1p= [1.146, 1.274]
list_QBDT3p= [1.110, 1.315]
list_QBDT7p= [1.088, 1.195]


list_QBDT = [list_QBDT1p, list_QBDT3p, list_QBDT7p]


list_dmu_GBDT = []
list_dmu_QBDT = []

for i in range(1, len(list_GBDT)):
   dmu = math.sqrt(pow(list_GBDT[i],2) - pow(list_GBDT[0],2))
   #list_dmu_GBDT.append(dmu*dmu)
   list_dmu_GBDT.append(dmu)
   dmu = math.sqrt(pow(list_QBDT[i-1][1],2) - pow(list_QBDT[i-1][0],2))
   #list_dmu_QBDT.append(dmu*dmu)
   list_dmu_QBDT.append(dmu)


for i in range(3):
   print 'dmu2', i , ':', list_dmu_GBDT[i], list_dmu_QBDT[i], list_dmu_QBDT[i]/list_dmu_GBDT[i]
