# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 00:30:10 2021
ver=1.0  
@author: NS
"""

import pandas as pd
import numpy as np
from gurobipy import *
from decimal import getcontext
import sys
from itertools import product
import plotly.graph_objects as go
import json
from os import path

getcontext().prec = 4
#distance read from excel
file='flcd.xls'
xl=pd.ExcelFile(file)
dist=xl.parse('Fixed_link_cost')
dist=np.array(dist)

#demand read from excel
file='hfa.xls'
xl=pd.ExcelFile(file)
dem=xl.parse('Fixed_link_cost')
dem=np.array(dem)


rho=float(sys.argv[1])
alpha= float(sys.argv[2])
perc=float(sys.argv[3])
fixcos=float(sys.argv[4])
capl=float(sys.argv[5])
capa=float(sys.argv[6])

nnode=len(dem[1,:])
#nnode=13
Cset = range(0,nnode)


o={}
c=dist/25000
ce={}
say1=0
say2=0
say3=0
for i in Cset:
    for k in Cset:
        if dist[i,k]<1000:
            ce[i,k]=dist[i,k]*1.39/1000000
            say1+=1
        elif dist[i,k]<3700:
            ce[i,k]=dist[i,k]*0.71/1000000
            say2+=1
        else:
           ce[i,k]=dist[i,k]*0.56/1000000
           say3+=1

ct={}
for i in Cset:
    for j in Cset:
        for k in Cset:
            ct[i,j,k]=ce[i,j]+ce[j,k]
    
for i in range(nnode):
    o[i]=0
    for k in range(nnode):
        o[i]=o[i]+dem[i,k]
demsum=sum(o[i] for i in Cset)        
d={}        
for k in range(nnode):
    d[k]=0
    for i in range(nnode):
        d[k]=d[k]+dem[i,k]
        
        

A={}
distmax=dist.max()
for i in Cset:
    for k in Cset:
        if dist[i,k]< perc*distmax:
            A[i,k]=1
        else:
            A[i,k]=0



fixedCost={} 
for i in Cset:
                fixedCost[i]=fixcos
detmas= Model("detmas")
detmas.params.outputflag=0

#decision variables
hopen = {} 
for p in Cset:
    hopen[p] = detmas.addVar(vtype=GRB.BINARY, obj=fixedCost[p], name="x%s" % p)
z={}
for i in Cset:
    for k in Cset:
        z[i,k]=detmas.addVar(obj=c[i,k] ,name='z%s.%s' % (i,k))
y={}
for i in Cset:
    for k in Cset:
        for l in Cset:
            y[i,k,l]=detmas.addVar(obj=alpha*c[k,l],name='y%s.%s.%s' % (i,k,l))
q={}
for i in Cset:
    for l in Cset:
        for j in Cset:
            q[i,l,j]=detmas.addVar(obj=c[l,j],name='q%s.%s.%s' % (i,l,j))

                
#constraints
sumdem={}
coll={}
for i in Cset:
    coll[i]=sum(z[i,k] for k in Cset if A[i,k]==1)
    sumdem[i]=detmas.addConstr(coll[i] == quicksum(dem[i,j] for j in Cset), "sumdem%s" % i)
    
#q demand ensure
qsumdem={}
quell={}
for i in Cset:
    for j in Cset:
        quell[i,j]=sum(q[i,l,j] for l in Cset if A[i,l]*A[l,j]==1)
        qsumdem[i,j]=detmas.addConstr(quell[i,j] ==dem[i,j], "qdem%s.%s" % (i,j) )

#flow ensure constraint
floens={}
y1ell={}
y2ell={}
q2ell={}
zwell={}
for i in Cset:
    for k in Cset:
        y1ell[i,k]=sum(y[i,k,l] for l in Cset if A[i,k]*A[k,l]==1)
        y2ell[i,k]=sum(y[i,l,k] for l in Cset if A[i,l]*A[l,k]==1)
        q2ell[i,k]=sum(q[i,k,j] for j in Cset if A[i,k]*A[k,j]==1)
        if A[i,k]==1:
            floens[i,k]=detmas.addConstr(y1ell[i,k]+q2ell[i,k]-y2ell[i,k]-z[i,k]==0, "ensco%s.%s" % (i,k))
        
#coverage constraint1
covcon={}
for l in Cset:
    for j in Cset:
        if A[l,j]==1:
            covcon[l,j]=detmas.addConstr(quicksum(q[i,l,j] for i in Cset) <= quicksum(dem[i,j] for i in Cset)*hopen[l], "covc%s.%s" % (l,j))
        
#coverage constraint2
covcot={}
for i in Cset:
    for k in Cset:
        if A[i,k]==1:
            covcot[i,k]=detmas.addConstr(z[i,k]<=quicksum(dem[i,j]*hopen[k] for j in Cset), "covct%s.%s" % (i,k))

otop=sum(o[i] for i  in Cset)            
capconl1={}
for i in Cset:
    for k in Cset:
        if A[i,k]==1:
            capconl1[k]=detmas.addConstr(z[i,k] <= ((1-hopen[i])*capl+ otop*hopen[i]) )
        
capconl2={}
for l in Cset:
    for j in Cset:
        if A[l,j]==1:
            capconl2[l,j]=detmas.addConstr(quicksum(q[i,l,j] for i in Cset) <= ((1-hopen[j])*capl+ otop*hopen[j]))
         
capcon={}
for k in Cset:
        capcon[k]=detmas.addConstr(quicksum(z[i,k] for i in Cset if A[i,k]==1) <= capa*8540006)             
 
detmas.modelSense = GRB.MINIMIZE
detmas.Params.MIPGap=1e-6
detmas.update() 
detmas.optimize()

numhum={}  
i=0            
for p in Cset:
    if round(hopen[p].x)!=0:
        numhum[p]=round(hopen[p].x)

      

carbrel={}
carbem=0
for i in Cset:
    carbrel[i]=sum(ce[i,r]*z[i,r].x for r in Cset)+sum(ct[i,sa1,sa2]*(y[i,sa1,sa2].x+q[i,sa1,sa2].x) for sa1,sa2 in product(Cset, repeat=2) )
    for k in Cset:
        carbem+=ce[i,k]*z[i,k].x
        for l in Cset:
          carbem+=ct[i,k,l]*(y[i,k,l].x+q[i,k,l].x)  
          
print("Writing results for instance {}")

new_data = {}
new_data["MaximumLinkCapacity"] = capl
new_data["MaximumHubCapacityRate"] = capa
new_data["CO2EmissionUnitCost"] = rho
new_data["InterhubDiscountFactor"] = alpha
new_data["CoveringRadius"] = perc
new_data["HubOpeningCost"] = fixcos
new_data["TotalCost"] = detmas.objVal
new_data["TotalWithoutCO2Cost"] = detmas.objVal-rho*carbem
new_data["TotalCarbonEmission"] = carbem
new_data["OpenedHubLocs"] = numhum
new_data["DemandCo2Emission"] = carbrel
new_data["RunTime"] = detmas.runtime
fname = "det_carb_hlco.json"

if path.exists(fname):
    with open(fname, 'r' ) as f:
        data = json.load(f)
        data.append(new_data)
else:
    data = [new_data]

with open(fname, 'w') as outfile:
    json.dump(data, outfile)