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
import os
import csv


def write_to_file(fname, data_dict, is_header):
    with open(fname, 'a') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=list(data_dict.keys()))
        if is_header == 1:
            writer.writeheader()
        writer.writerow(data_dict)

getcontext().prec = 4
# distance read from excel
file = 'flcd.xls'
xl = pd.ExcelFile(file)
dist = xl.parse('Fixed_link_cost')
dist = np.array(dist)

# demand read from excel
file = 'hfa.xls'
xl = pd.ExcelFile(file)
dem = xl.parse('Fixed_link_cost')
dem = np.array(dem)




nnode = len(dem[1, :])
# nnode=13
Cset = range(0, nnode)


o = {}
c = dist/25000
ce = {}
say1 = 0
say2 = 0
say3 = 0
for i in Cset:
    for k in Cset:
        if dist[i, k] < 1000:
            ce[i, k] = dist[i, k]*1.39/1000000
            say1 += 1
        elif dist[i, k] < 3700:
            ce[i, k] = dist[i, k]*0.71/1000000
            say2 += 1
        else:
            ce[i, k] = dist[i, k]*0.56/1000000
            say3 += 1

ct = {}
for i in Cset:
    for j in Cset:
        for k in Cset:
            ct[i, j, k] = ce[i, j]+ce[j, k]

for i in range(nnode):
    o[i] = 0
    for k in range(nnode):
        o[i] = o[i]+dem[i, k]
demsum = sum(o[i] for i in Cset)
d = {}
for k in range(nnode):
    d[k] = 0
    for i in range(nnode):
        d[k] = d[k]+dem[i, k]

params = np.loadtxt("carb_paramshlc.csv", delimiter=",", skiprows=0)
all_results=[]
result_file = os.path.join("results", "result_hlc_w.csv")        

for case_idx, param in enumerate(params):
    rho, alpha, perc, fixcos, capl, capa = param 
    A = {}
    distmax = dist.max()
    for i in Cset:
        for k in Cset:
            if dist[i, k] < perc*distmax:
                A[i, k] = 1
            else:
                A[i, k] = 0
    
    
    fixedCost = {}
    for i in Cset:
        fixedCost[i] = fixcos
    detmas = Model("detmas")
    detmas.params.outputflag = 0
    
    # decision variables
    hopen = {}
    for p in Cset:
        hopen[p] = detmas.addVar(
            vtype=GRB.BINARY, obj=fixedCost[p], name="x%s" % p)
    z = {}
    for i in Cset:
        for k in Cset:
            z[i, k] = detmas.addVar(
                obj=c[i, k]+rho*ce[i, k], name='z%s.%s' % (i, k))
    y = {}
    for i in Cset:
        for k in Cset:
            for l in Cset:
                y[i, k, l] = detmas.addVar(
                    obj=alpha*c[k, l]+rho*ct[i, k, l], name='y%s.%s.%s' % (i, k, l))
    q = {}
    for i in Cset:
        for l in Cset:
            for j in Cset:
                q[i, l, j] = detmas.addVar(
                    obj=c[l, j]+rho*ct[i, l, j], name='q%s.%s.%s' % (i, l, j))
    
    
    # constraints
    sumdem = {}
    coll = {}
    for i in Cset:
        coll[i] = sum(z[i, k] for k in Cset if A[i, k] == 1)
        sumdem[i] = detmas.addConstr(coll[i] == quicksum(
            dem[i, j] for j in Cset), "sumdem%s" % i)
    
    # q demand ensure
    qsumdem = {}
    quell = {}
    for i in Cset:
        for j in Cset:
            quell[i, j] = sum(q[i, l, j] for l in Cset if A[i, l]*A[l, j] == 1)
            qsumdem[i, j] = detmas.addConstr(
                quell[i, j] == dem[i, j], "qdem%s.%s" % (i, j))
    
    # flow ensure constraint
    floens = {}
    y1ell = {}
    y2ell = {}
    q2ell = {}
    zwell = {}
    for i in Cset:
        for k in Cset:
            y1ell[i, k] = sum(y[i, k, l] for l in Cset if A[i, k]*A[k, l] == 1)
            y2ell[i, k] = sum(y[i, l, k] for l in Cset if A[i, l]*A[l, k] == 1)
            q2ell[i, k] = sum(q[i, k, j] for j in Cset if A[i, k]*A[k, j] == 1)
            if A[i, k] == 1:
                floens[i, k] = detmas.addConstr(
                    y1ell[i, k]+q2ell[i, k]-y2ell[i, k]-z[i, k] == 0, "ensco%s.%s" % (i, k))
    
    # coverage constraint1
    covcon = {}
    for l in Cset:
        for j in Cset:
            if A[l, j] == 1:
                covcon[l, j] = detmas.addConstr(quicksum(q[i, l, j] for i in Cset) <= quicksum(
                    dem[i, j] for i in Cset)*hopen[l], "covc%s.%s" % (l, j))
    
    # coverage constraint2
    covcot = {}
    for i in Cset:
        for k in Cset:
            if A[i, k] == 1:
                covcot[i, k] = detmas.addConstr(z[i, k] <= quicksum(
                    dem[i, j]*hopen[k] for j in Cset), "covct%s.%s" % (i, k))
    
    otop = sum(o[i] for i in Cset)
    capconl1 = {}
    for i in Cset:
        for k in Cset:
            if A[i, k] == 1:
                capconl1[k] = detmas.addConstr(
                    z[i, k] <= ((1-hopen[i])*capl + otop*hopen[i]))
    
    capconl2 = {}
    for l in Cset:
        for j in Cset:
            if A[l, j] == 1:
                capconl2[l, j] = detmas.addConstr(
                    quicksum(q[i, l, j] for i in Cset) <= ((1-hopen[j])*capl + otop*hopen[j]))
    
    capcon = {}
    for k in Cset:
        capcon[k] = detmas.addConstr(quicksum(z[i, k]
                                              for i in Cset if A[i, k] == 1) <= capa*8540006)
    
    detmas.modelSense = GRB.MINIMIZE
    detmas.Params.MIPGap = 1e-6
    detmas.update()
    detmas.optimize()
    
    numhum = {}
    i = 0
    for p in Cset:
        if round(hopen[p].x) != 0:
            numhum[p] = round(hopen[p].x)
    
    
    carbrel = {}
    carbem = 0
    for i in Cset:
        carbrel[i] = sum(ce[i, r]*z[i, r].x for r in Cset)+sum(ct[i, sa1, sa2] *
                                                               (y[i, sa1, sa2].x+q[i, sa1, sa2].x) for sa1, sa2 in product(Cset, repeat=2))
        for k in Cset:
            carbem += ce[i, k]*z[i, k].x
            for l in Cset:
                carbem += ct[i, k, l]*(y[i, k, l].x+q[i, k, l].x)

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
    all_results.append(new_data)
    write_to_file(result_file, new_data, case_idx==0) 
    
    # print("Writing results for instance {}")
    
    # df = pd.read_csv(r"C:\Users\Administrator\Desktop\Python\locations.csv")
    # opened_hub_nos = list(numhum.keys())   
    # #print(opened_hub_nos)
    # df_draw_squares = df.loc[opened_hub_nos, :]   # bunlar yazdrlacak locasyonlar
    
    
    
    # df['text'] = df['airport'] + '' + df['city'] + ', ' + df['state'] + '' + 'Arrivals: ' 
    # df['cnt']=df.index.to_series().map(carbrel)
    # df['cnt'].astype(str)
    # fig = go.Figure()
    
    # fig.add_trace(go.Scattergeo(
    #         locationmode = 'USA-states',
    #         lon = df['long'],
    #         lat = df['lat'],
    #         showlegend=False,
    #         text = df['text'],
    #         mode = 'markers',
    #         marker = dict(
    #             size = 8,
    #             opacity = 0.8,
    #             reversescale = True,
    #             autocolorscale = False,
    #             symbol = 'square',
    #             line = dict(
    #                 width=0,
    #                 color='rgba(102, 102, 102)'
    #             ),
    #             colorscale = [[0, 'rgb(0,60,0)'], [0.5, 'rgb(0,128,0)'], [1, 'rgb(0,255,0)']],
    #             cmin = 0,
    #             color = df['cnt'],
    #             cmax = 2300 ,
    #             colorbar_title="" )))
    
    # df_draw_squares['cnt'] = [0]*len(df_draw_squares)
    # print(df_draw_squares)
    # fig.add_trace(go.Scattergeo(
    #         name="Opened Hubs",
    #         locationmode = 'USA-states',
    #         lon = df_draw_squares['long'],
    #         lat = df_draw_squares['lat'],
    #         mode = 'markers',
    #         marker = dict(
    #             size = 11.0,
    #             opacity = 1.0,
    #             reversescale = False,
    #             autocolorscale = False,
    #             color='green',
    #             symbol = 'circle-open',    
              
    # print("Writing results for instance {}")
    #             line = dict(
    #                 width=1,
    #                 color='green'
    #             ))))
    
    
    # fig.update_layout(
    #         geo = dict(
    #             scope='usa',
    #             projection_type='albers usa',
    #             showland = True,
    #             landcolor = "rgb(250, 250, 250)",
    #             subunitcolor = "rgb(217, 217, 217)",
    #             countrycolor = "rgb(217, 217, 217)",
    #             countrywidth = 0.5,
    #             subunitwidth = 0.5
    #         ),  annotations=[dict(
    #           # Don't specify y position,because yanchor="middle" do it for you
    #           x=1.02,
    #           align="right",
    #           valign="top",
    #           text='Carbon Emissions (kg/year)',
    #           showarrow=False,
    #           xref="paper",
    #           yref="paper",
    #           xanchor="right",
    #           yanchor="middle",
    #           # Parameter textangle allow you to rotate annotation how you want
    #           textangle=-90
    #         )
    #     ], 
    #         legend=dict(
    #                 yanchor="top",
    #                 y=0.99,
    #                 xanchor="center",
    #                 x=0.50)
    #     )
    # fig.write_image('images/hlc/'+str(capl)+str(capa)+str(rho) +
    #                 str(alpha)+str(perc)+str(fixcos)+'.pdf')

result_file = os.path.join('results', 'carb_hlc_w.json')
with open(result_file, 'w') as outfile:
    json.dump(all_results, outfile)   