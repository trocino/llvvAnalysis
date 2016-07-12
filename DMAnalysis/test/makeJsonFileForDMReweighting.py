#!/usr/bin/env python
import os,sys
import json
import getopt
import commands
from math import sqrt, log10

print '{'
print '  "proc": ['

mds = ['V', 'A'] 

gqs = ['', '_GQ0p25'] 

mxs = [1, 5, 10, 25, 50, 60, 70, 80, 90, 100, 110, 125, 135, 150, 175, 200, 250, 300, 400, 500, 600, 700, 800, 900, 1000]

mvs = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 300, 325, 350, 400, 450, 525, 600, 725, 800, 925, 1000, 1125, 1200, 1325, 1400, 1525, 1600, 1725, 1800, 1925, 2000, 2500, 3000, 3500, 4000, 5000]


vpts = [
    ["MC13TeV_DM_V_Mx1Mv10"     ,    1,   10, 22.926	       ],
    ["MC13TeV_DM_V_Mx1Mv20"     ,    1,   20,  9.1442	       ],
    ["MC13TeV_DM_V_Mx1Mv50"     ,    1,   50,  2.4785	       ],
    ["MC13TeV_DM_V_Mx1Mv100"    ,    1,  100,  0.81966	       ],
    ["MC13TeV_DM_V_Mx1Mv200"    ,    1,  200,  0.24152	       ],
    ["MC13TeV_DM_V_Mx1Mv300"    ,    1,  300,  0.107078	       ],
    ["MC13TeV_DM_V_Mx1Mv500"    ,    1,  500,  0.027291	       ],
    ["MC13TeV_DM_V_Mx1Mv1000"   ,    1, 1000,  0.0037816       ],
    ["MC13TeV_DM_V_Mx1Mv2000"   ,    1, 2000,  0.00030459      ],
    ["MC13TeV_DM_V_Mx1Mv5000"   ,    1, 5000,  0.0000049242    ],
    ["MC13TeV_DM_V_Mx10Mv10"    ,   10,   10,  1.94108	       ],
    ["MC13TeV_DM_V_Mx10Mv20"    ,   10,   20,  4.6803	       ],
#    ["MC13TeV_DM_V_Mx10Mv50"    ,   10,   50,  2.4185	       ],
    ["MC13TeV_DM_V_Mx10Mv100"   ,   10,  100,  0.81682	       ],
    ["MC13TeV_DM_V_Mx10Mv5000"  ,   10, 5000,  0.0000049257    ],
    ["MC13TeV_DM_V_Mx50Mv10"    ,   50,   10,  0.114921	       ],
    ["MC13TeV_DM_V_Mx50Mv50"    ,   50,   50,  0.144129	       ],
    ["MC13TeV_DM_V_Mx50Mv95"    ,   50,   95,  0.341923	       ],
    ["MC13TeV_DM_V_Mx50Mv200"   ,   50,  200,  0.22546	       ],
    ["MC13TeV_DM_V_Mx50Mv300"   ,   50,  300,  0.104327	       ],
    ["MC13TeV_DM_V_Mx50Mv5000"  ,   50, 5000,  0.0000049063    ],
    ["MC13TeV_DM_V_Mx150Mv10"   ,  150,   10,  0.0112056       ],
    ["MC13TeV_DM_V_Mx150Mv200"  ,  150,  200,  0.0189941       ],
    ["MC13TeV_DM_V_Mx150Mv295"  ,  150,  295,  0.04168	       ],
    ["MC13TeV_DM_V_Mx150Mv500"  ,  150,  500,  0.022839	       ],
    ["MC13TeV_DM_V_Mx150Mv5000" ,  150, 5000,  0.0000046839    ],
    ["MC13TeV_DM_V_Mx500Mv10"   ,  500,   10,  0.00027615      ],
#    ["MC13TeV_DM_V_Mx500Mv500"  ,  500,  500,  0.000383387     ],
    ["MC13TeV_DM_V_Mx500Mv995"  ,  500,  995,  0.0010751174    ],
    ["MC13TeV_DM_V_Mx500Mv2000" ,  500, 2000,  0.0002141       ],
    ["MC13TeV_DM_V_Mx500Mv5000" ,  500, 5000,  0.0000031352    ],
    ["MC13TeV_DM_V_Mx1000Mv10"  , 1000,   10,  0.0000109818    ],
    ["MC13TeV_DM_V_Mx1000Mv1000", 1000, 1000,  0.00001590609   ],
    ["MC13TeV_DM_V_Mx1000Mv1995", 1000, 1995,  0.00004794683192],
    ["MC13TeV_DM_V_Mx1000Mv5000", 1000, 5000,  0.00000143883   ]
    ]

apts = [
    ["MC13TeV_DM_A_Mx1Mv10"	,    1,   10, 24.084           ], 
    ["MC13TeV_DM_A_Mx1Mv20"	,    1,   20,  9.5741          ], 
#    ["MC13TeV_DM_A_Mx1Mv50"	,    1,   50,  2.4991          ], 
#    ["MC13TeV_DM_A_Mx1Mv100"    ,    1,  100,  0.82149	       ], 
    ["MC13TeV_DM_A_Mx1Mv200"    ,    1,  200,  0.24167	       ], 
    ["MC13TeV_DM_A_Mx1Mv300"    ,    1,  300,  0.107074	       ], 
    ["MC13TeV_DM_A_Mx1Mv500"    ,    1,  500,  0.030356	       ], 
#    ["MC13TeV_DM_A_Mx1Mv1000"   ,    1, 1000,  0.0039114       ], 
    ["MC13TeV_DM_A_Mx1Mv2000"   ,    1, 2000,  0.00030726      ], 
    ["MC13TeV_DM_A_Mx1Mv5000"   ,    1, 5000,  0.0000049292    ], 
    ["MC13TeV_DM_A_Mx10Mv10"    ,   10,   10,  1.02349	       ], 
    ["MC13TeV_DM_A_Mx10Mv20"    ,   10,   20,  1.85535	       ], 
    ["MC13TeV_DM_A_Mx10Mv50"    ,   10,   50,  1.96865	       ], 
    ["MC13TeV_DM_A_Mx10Mv100"   ,   10,  100,  0.77058	       ], 
    ["MC13TeV_DM_A_Mx10Mv5000"  ,   10, 5000,  0.0000049256    ], 
    ["MC13TeV_DM_A_Mx50Mv10"    ,   50,   10,  0.057948	       ], 
    ["MC13TeV_DM_A_Mx50Mv50"    ,   50,   50,  0.0673412       ], 
    ["MC13TeV_DM_A_Mx50Mv95"    ,   50,   95,  0.1206166       ], 
    ["MC13TeV_DM_A_Mx50Mv200"   ,   50,  200,  0.1581	       ], 
    ["MC13TeV_DM_A_Mx50Mv300"   ,   50,  300,  0.08755	       ], 
#    ["MC13TeV_DM_A_Mx50Mv5000"  ,   50, 5000,  0.0000047845    ], 
    ["MC13TeV_DM_A_Mx150Mv10"   ,  150,   10,  0.0048086       ], 
    ["MC13TeV_DM_A_Mx150Mv200"  ,  150,  200,  0.0069606       ], 
    ["MC13TeV_DM_A_Mx150Mv295"  ,  150,  295,  0.01262797      ], 
    ["MC13TeV_DM_A_Mx150Mv500"  ,  150,  500,  0.0148225       ], 
    ["MC13TeV_DM_A_Mx150Mv5000" ,  150, 5000,  0.000004152     ], 
    ["MC13TeV_DM_A_Mx500Mv10"   ,  500,   10,  0.000089041     ], 
    ["MC13TeV_DM_A_Mx500Mv500"  ,  500,  500,  0.000114785     ], 
    ["MC13TeV_DM_A_Mx500Mv995"  ,  500,  995,  0.00028419647   ], 
    ["MC13TeV_DM_A_Mx500Mv2000" ,  500, 2000,  0.000125743     ], 
    ["MC13TeV_DM_A_Mx500Mv5000" ,  500, 5000,  0.0000020417    ], 
    ["MC13TeV_DM_A_Mx1000Mv10"  , 1000,   10,  0.0000027221    ], 
    ["MC13TeV_DM_A_Mx1000Mv1000", 1000, 1000,  0.00000368319   ], 
    ["MC13TeV_DM_A_Mx1000Mv1995", 1000, 1995,  0.00001021123579], 
    ["MC13TeV_DM_A_Mx1000Mv5000", 1000, 5000,  6.863E-7	       ]
    ]



cnt = 0

for imd in mds: 
    for igq in gqs: 
        for imx in mxs:
            for imv in mvs:
                idtag = 'MC13TeV_DM_'+imd+igq+'_Mx'+str(imx)+'Mv'+str(imv) 
                #if idtag in vpts or idtag in apts: continue 
                isfullsim = False 
                mindiff = 999999. 
                minipt = ["", -1, -1, -1.] 
                xpts = [] 
                if   imd == 'V': xpts = vpts 
                elif imd == 'A': xpts = apts 
                for ipt in xpts: 
                    if idtag == ipt[0]: 
                        isfullsim = True 
                        break 
                    elif imx > ipt[1] or imv > ipt[2]: 
                        continue 
                    else: 
                        diff = sqrt( (log10(ipt[1]/imx))**2 + (log10(ipt[2]/imv))**2 ) 
                        if diff < mindiff: 
                            mindiff = diff 
                            minipt = ipt 
                if isfullsim: continue 
                if minipt[3] < 0.0: 
                    print " *** WARNING: minimization failed for point (mx, MV) = (" + str(imx) + ", " + str(imv) + ") ***"  
                    continue 

                if cnt == 0:
                    cnt = 1 
                else: 
                    print '    },'
                print '    {' 
                print '      "isdata": false,'
                print '      "color": "809",'
                print '      "issignal": true,'
                print '      "tag": "DM('+str(imx)+')M'+imd+'('+str(imv)+')",'
                print '      "lwidth": 2,' 
                print '      "marker": 809,'
                print '      "data": ['
                print '        {'
                print '          "dtag": "'+idtag+'",'
                print '          "xsec": '+str(minipt[3])+','
                print '          "split": 1,'
                print '          "br": ['
                print '            1.0'
                print '          ]'
                print '        }'
                print '      ],'
                print '      "fill": 0,'
                print '      "spimpose": true,'
                print '      "lcolor": "809",'
                print '      "lstyle": 8,'
                print '      "doWIMPreweighting": true'

print '    }'
print '  ]'
print '}'

exit 
