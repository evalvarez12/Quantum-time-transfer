# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')
plt.rcParams["font.family"] = "Times New Roman"

##################################
## OLD DATA

# ss20 = [168080.02362373075,
#  139312.93103616606,
#  110397.11912910634,
#  87582.82707938358,
#  68483.83089332229,
#  42188.47880154584,
#  16141.194982485085,
#  7603.808323072739,
#  2976.236762023574,
#  0.0,
#  0.0]

# qs20 =[0.06296370472301642,
#  0.06354876358578446,
#  0.06475344663709771,
#  0.06645050053308645,
#  0.06786098210062304,
#  0.0708113167613888,
#  0.07668211943111482,
#  0.08939633845614606,
#  0.12360414813612813,
#  0.2115326272480434,
#  0.7333504460787376]


# ss_ES20 = [416933.02074454713,
#  337732.02074454713,
#  267614.02074454713,
#  209319.0207445471,
#  160996.0207445471,
#  96513.02074454712,
#  56784.02074454712,
#  26795.020744547117,
#  10193.020744547115,
#  4075.0207445471156,
#  1688.0207445471158]

# qs_ES20 = [0.051177760270251955,
#  0.051312864335332706,
#  0.05149477651569437,
#  0.0517022685419102,
#  0.051964819840179526,
#  0.05262070138707676,
#  0.053520703915799266,
#  0.055464537241273135,
#  0.06012205602188101,
#  0.06927730476651459,
#  0.08515919477228262]


# ss100 = [783757.5449182872,
#  641841.5626828604,
#  513629.71988776035,
#  406746.39956781577,
#  315972.2426138987,
#  192655.93688537445,
#  117452.0619768067,
#  58286.611810268114,
#  14951.370177194656,
#  0.0,
#  0.0]


# qs100 = [0.05921149040696619,
#  0.059688541983419996,
#  0.06006485931668837,
#  0.060533458858774876,
#  0.061250313390559266,
#  0.0629952207544983,
#  0.06527063564304592,
#  0.06971236017010961,
#  0.08297948152125371,
#  0.11043397725705788,
#  0.1859861106778129]

# ss400 = [2998531.46987778,
#  2458015.729328914,
#  1954819.6078345277,
#  1542660.6742851653,
#  1201876.3020162424,
#  726548.1615849381,
#  441214.7166826235,
#  217735.24157367583,
#  91807.75003566957,
#  32720.324643622575,
#  0.0]

# qs400 = [0.057471618122206775,
#####################################


dists = [7000, 8000, 9000, 10000, 11000, 13000, 15000, 18000, 22000, 26000, 30000]


################ Nighttime 

ss20_l = [197700, 162363, 129118, 101957, 80511, 48195,
          29234, 7302, 2817]

ss20_u = [423748, 348495, 279291, 222108, 103267, 62758,
          38051, 8400, 6048]

qs20_l = [0.0631, 0.0638, 0.0646, 0.0656, 0.0672, 0.0715,
          0.0767, 0.0876, 0.1249]

qs20_u = [0.0951, 0.0876, 0.0931, 0.0916, 0.0993, 0.0963,
          0.1028, 0.1130, 0.1301]

ss_e = [416336, 340391, 267397, 210614, 161914, 95789,
        56911, 26828, 10031, 4069.0207445471156, 1777.0207445471158]

qs_e = [0.0511, 0.0513, 0.0514, 0.0516, 0.0519, 0.0526,
        0.0535, 0.0554, 0.0601, 0.06892400699932003, 0.08533041812611217]

ss400_l = [3516922, 2888315.0, 2299482, 1818229,
           1415230, 860970, 523590, 256304,
           107194, 29642, 15954]

ss400_u = [7774613, 9245584, 7355002, 5350911, 
           3767015, 2358578, 1364198, 649559, 
           203513, 67035, 21216]

qs400_l = [0.0574, 0.0578, 0.0580, 0.0584,
           0.0588, 0.0596, 0.0609, 0.0637, 
           0.0709, 0.0870, 0.1230]

qs400_u = [0.0979, 0.0947, 0.0945, 0.0998, 
           0.0956, 0.1017, 0.1002, 0.0908, 
           0.0889, 0.1145, 0.1374]


################ Daytime 


ss20_l_d = [40775, 29484]

ss20_u_d = [63881, 30001]

qs20_l_d = [0.1086, 0.1352]

qs20_u_d = [0.1279, 0.1359]

ss_e_d = [86582.03722735579, 56338.03722735579, 37137.03722735579, 25744.037227355788,
          18441.037227355788, 10900.037227355788, 7255.037227355788, 4977.037227355788,
          3904.0372273557873, 3492.0372273557873, 3327.0372273557873]

qs_e_d = [0.06916051599756254, 0.07882349971087002, 0.09302295521572022,
          0.11136727795086988, 0.13495587231653183, 0.192197070937214,
          0.2620681944263564, 0.3569110571930854, 0.43913261760012434,
          0.48378631557652396, 0.5046814688571426]

ss400_l_d = [904996, 656263]

ss400_u_d = [1520276, 862188]

qs400_l_d = [0.0921, 0.1111]

qs400_u_d = [0.1177, 0.1237]



########### confidence

ss90conf = [193262.0, 159208.0, 126531.0, 100013.0, 77851.0, 47555.0,
            28914.999999999996, 7371.999999999999, 2826.9999999999995, 1384.0, 0.0]
qs90conf = [0.061321875226214254, 0.06231669528232329, 0.0626981536556063,
            0.06406294149916791, 0.06508827467310813, 0.06785031100557534,
            0.07213492671474002, 0.08035868886955023, 0.1032166112019585,
            0.16341887222037807, 0.35940611700170294]

ss50conf = [188518.0, 153884.0, 122790.0, 97322.0, 76471.0, 46255.00000000001,
            27607.0, 13863.000000000002, 2871.0, 1432.0, 0.0]
qs50conf = [0.05979842746296182, 0.06062443158199925, 0.06120862163075198,
            0.06153993484642844, 0.06253761941561266, 0.06448754679780992,
            0.0676360948079931, 0.07388173900951583, 0.09168600981096839,
            0.12450709381238019, 0.23794013811594503]

#################################################################
### PLOTS

cs = plt.rcParams['axes.prop_cycle'].by_key()['color']



fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)


ax.plot(dists[:-2], ss20_l, color=cs[0])
ax.plot(dists[:-2], ss20_u, color=cs[0])
ax.fill_between(dists[:-2], ss20_l, ss20_u, alpha=0.5, color=cs[0], label='20 MHz')

ax.plot(dists[:], ss400_l, color=cs[1])
ax.plot(dists[:], ss400_u, color=cs[1])
ax.fill_between(dists[:], ss400_l, ss400_u, alpha=0.5, color=cs[1], label='400 MHz')


ax.plot(dists, ss_e, color=cs[2], label='ES - 20 Mhz')

ax.set_yscale('log')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Coincidences')

plt.legend()


fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)


ax.plot(dists[:-2], qs20_l, color=cs[0])
ax.plot(dists[:-2], qs20_u, color=cs[0])
ax.fill_between(dists[:-2], qs20_l, qs20_u, alpha=0.5, color=cs[0], label='20 MHz')

ax.plot(dists[:], qs400_l, color=cs[1])
ax.plot(dists[:], qs400_u, color=cs[1])
ax.fill_between(dists[:], qs400_l, qs400_u, alpha=0.5, color=cs[1], label='400 MHz')


ax.plot(dists, qs_e, color=cs[2], label='ES - 20 Mhz')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('QBER')


plt.legend()



##### Daytime plot

fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)


ax.plot(dists[:2], ss20_l_d, color=cs[0])
ax.plot(dists[:2], ss20_u_d, color=cs[0])
ax.fill_between(dists[:2], ss20_l_d, ss20_u_d, alpha=0.5, color=cs[0], label='20 MHz')

ax.plot(dists[:2], ss400_l_d, color=cs[1])
ax.plot(dists[:2], ss400_u_d, color=cs[1])
ax.fill_between(dists[:2], ss400_l_d, ss400_u_d, alpha=0.5, color=cs[1], label='400 MHz')


ax.plot(dists[:5], ss_e_d[:5], color=cs[2], label='ES - 20 Mhz')

ax.set_yscale('log')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Coincidences')

plt.legend()


fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)


ax.plot(dists[:2], qs20_l_d, color=cs[0])
ax.plot(dists[:2], qs20_u_d, color=cs[0])
ax.fill_between(dists[:2], qs20_l_d, qs20_u_d, alpha=0.5, color=cs[0], label='20 MHz')

ax.plot(dists[:2], qs400_l_d, color=cs[1])
ax.plot(dists[:2], qs400_u_d, color=cs[1])
ax.fill_between(dists[:2], qs400_l_d, qs400_u_d, alpha=0.5, color=cs[1], label='400 MHz')


ax.plot(dists[:5], qs_e_d[:5], color=cs[2], label='ES - 20 Mhz')

ax.set_xlabel('Distance (m)')
ax.set_ylabel('QBER')
plt.legend()




########### Confidence plot
fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)


ax.plot(dists[:-2], ss20_l, color=cs[0])
ax.plot(dists[:-1], ss90conf[:-1], color=cs[1])
ax.plot(dists[:-1], ss50conf[:-1], color=cs[2])



ax.set_yscale('log')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Coincidences')

plt.legend()


fig, ax = plt.subplots()
fig.set_size_inches(18.5*.3, 10.5*.3)


ax.plot(dists[:-2], qs20_l, color=cs[0], label='99% confidence')
ax.plot(dists[:-1], qs90conf[:-1], color=cs[1], label='90% confidence')
ax.plot(dists[:-1], qs50conf[:-1], color=cs[2], label='50% confidence')


ax.set_xlabel('Distance (m)')
ax.set_ylabel('QBER')


plt.legend()








plt.show()


# fig = plt.figure()
# fig.set_size_inches(18.5*.3, 10.5*.3)

# plt.plot(dists, ss20, linewidth=3, label='BB84 - 20 MHz')
# plt.plot(dists, ss100, linewidth=3, label='BB84 - 100 MHz')
# plt.plot(dists, ss400, linewidth=3, label='BB84 - 400 MHz')
# plt.plot(dists, ss_ES20, linewidth=3, label='ES - 20 MHz')


# plt.xlabel('Distance (m)',fontsize=13)
# plt.ylabel('Coincidences', fontsize=13)

# plt.yticks(fontsize=13)
# plt.xticks(fontsize=13)


# plt.legend()
# # fig.savefig('loss.pdf', dpi=200)

# fig = plt.figure()
# fig.set_size_inches(18.5*.3, 10.5*.3)

# plt.plot(dists[:-2], qs20[:-2], linewidth=3, label='BB84 - 20 MHz')
# plt.plot(dists[:-2], qs100[:-2], linewidth=3, label='BB84 - 100 MHz')
# plt.plot(dists[:-1], qs400[:-1], linewidth=3, label='BB84 - 400 MHz')
# plt.plot(dists, qs_ES20, linewidth=3, label='ES - 20 MHz')


# plt.xlabel('Distance (m)',fontsize=13)
# plt.ylabel('QBER', fontsize=13)

# plt.yticks(fontsize=13)
# plt.xticks(fontsize=13)


# plt.legend()


# ####################### 20 Mhz source rate - night
# ---> d =  7000
# x:  [0.55551414 0.099999   0.19923096]
#       -        
# Nm:  29103.0
# Nv:  21120.0
# Y1:  0.011610031239379833
# m safe:  True
# v safe:  True
# qs:  [0.11232985484776985, 0.07294666843087315, 0.12562881487103664, 0.03845151238315885]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  83964.02074454712
# QBER:  0.052820171751474444

# ---> d =  8000
# x:  [0.4 0.1 0.4]
#       -        
# Nm:  26540.0
# Nv:  9970.999999999998
# Y1:  0.006997336464345395
# m safe:  True
# v safe:  True
# qs:  [0.10143542431420519, 0.08149283548755716, 0.11896503105211469, 0.03131699736256704]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  53006.02074454712
# QBER:  0.05366515375754354

# ---> d =  9000
# x:  [0.4 0.1 0.4]
#       -        
# Nm:  17378.0
# Nv:  6288.000000000001
# Y1:  0.004109777732871815
# m safe:  True
# v safe:  True
# qs:  [0.1168290177738811, 0.09254542912857816, 0.1399518998315338, 0.02433748954327025]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  33982.02074454712
# QBER:  0.054754804755875144

# ---> d =  10000
# x:  [0.4 0.1 0.4]
#       -        
# Nm:  11554.0
# Nv:  4323.0
# Y1:  0.002743717540806905
# m safe:  True
# v safe:  True
# qs:  [0.12096144801117176, 0.10086655754781951, 0.1468027520229386, 0.017596231964104234]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  22589.020744547117
# QBER:  0.05608254019697867

# ---> d =  11000
# x:  [0.66825592 0.099999   0.13167172]
#       -        
# Nm:  4205.0
# Nv:  4172.0
# N total:  8377.0
# Y1:  0.001801441159034726
# m safe:  False
# v safe:  True
# qs:  [0.17990937658636943, 0.1038525573024503, 0.20418964274549134, 0.017653437254263122]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  15340.020744547115
# QBER:  0.057738297962556834

# ---> d =  13000
# x:  [0.75231686 0.11102412 0.13566445]
#       -        
# Nm:  2517.0
# Nv:  2321.0
# N total:  4838.0
# Y1:  0.0007721050137112032
# m safe:  False
# v safe:  True
# qs:  [0.25920285676960364, 0.1398538420042708, 0.2983842608706651, -0.006296494991024021]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  7710.020744547116
# QBER:  0.06212287013860186

# ---> d =  15000
# x:  [0.67158907 0.13462315 0.16711235]
#       -        
# Nm:  7339.0
# Nv:  7315.000000000001
# N total:  14654.0
# Y1:  0.00045119620830013935
# m safe:  False
# v safe:  True
# qs:  [0.1921812410151696, 0.11040753392655012, 0.22842859994521936, 0.011355501935160356]
# ------ ES
# signal:  4232.020744547116
# QBER:  0.06839710441806472


# ####################### 100 Mhz source rate - night

# ---> d =  7000
# x:  [0.6, 0.1, 0.4]
#       -        
# Nm:  314852.0
# Nv:  78801.99999999999
# N total:  393654.0
# Y1:  0.012050291059814065
# m safe:  True
# v safe:  True
# qs:  [0.10351183792440437, 0.06521136970354444, 0.11328631900889856, 0.04486495141743915]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  415979.1037227356
# QBER:  0.051456191071911614

# ---> d =  8000
# x:  [0.6, 0.1, 0.4]
#       -        
# Nm:  200571.00000000003
# Nv:  50446.0
# N total:  251017.00000000003
# Y1:  0.007618873075771082
# m safe:  True
# v safe:  True
# qs:  [0.10551640106986042, 0.06737691790291377, 0.11585838120508743, 0.043464520258498396]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  266378.1037227356
# QBER:  0.05192856964347389


# ---> d =  9000
# x:  [0.4, 0.1, 0.4]
#       -        
# Nm:  86169.0
# Nv:  32142.999999999996
# N total:  118312.0
# Y1:  0.004756354855992746
# m safe:  True
# v safe:  True
# qs:  [0.0917389520392612, 0.07058120403146027, 0.10470398354456302, 0.03987882601805384]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  170662.10372273557
# QBER:  0.0525801313130576

# ---> d =  10000
# x:  [0.5, 0.1, 0.4]
#       -        
# Nm:  72147.0
# Nv:  22017.0
# N total:  94164.0
# Y1:  0.003216302133681871
# m safe:  True
# v safe:  True
# qs:  [0.10150078507501398, 0.07347910509328762, 0.11408053338952263, 0.038602043502470644]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  114033.10372273558
# QBER:  0.05339165022534882

# ---> d =  11000
# x:  [0.5, 0.1, 0.4]
#       -        
# Nm:  48099.0
# Nv:  14470.999999999998
# N total:  62570.0
# Y1:  0.0020276385747890397
# m safe:  True
# v safe:  True
# qs:  [0.1099119288704751, 0.07955390257704101, 0.12497359149232883, 0.03460361576120647]
# ---- P_fail
# (0.0, 0.0, 0.0)
# Key rate:  -9.58288677484898e-05

# ---> d =  13000
# x:  [0.4, 0.1, 0.4]
#       -        
# Nm:  19274.0
# Nv:  7422.0
# N total:  26696.0
# Y1:  0.0009574407800448958
# m safe:  True
# v safe:  True
# qs:  [0.11353266360984021, 0.09356995542185394, 0.13687844394746, 0.020149542259361105]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  38768.10372273558
# QBER:  0.057422794463172035

# ---> d =  15000
# x:  [0.67202143 0.13500392 0.16612484]
#       -        
# Nm:  7387.999999999999
# Nv:  7284.0
# N total:  14672.0
# Y1:  0.0004398852504196671
# m safe:  False
# v safe:  True
# qs:  [0.1992869802738972, 0.1124474280738991, 0.23723975383448132, 0.010365973040082976]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  21017.103722735577
# QBER:  0.06199707209019091

# ---> d =  18000
# x:  [0.78374544 0.16956526 0.18620061]
#       -        
# Nm:  4221.999999999999
# Nv:  4132.0
# N total:  8354.0
# Y1:  0.00016244482300767916
# m safe:  False
# v safe:  True
# qs:  [0.3017922452647431, 0.16281068973318355, 0.3751478590382635, -0.03726380879411057]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  9667.10372273558
# QBER:  0.07289737309775479


# ####################### 400 Mhz source rate - night
# ---> d =  7000
# x:  [0.6, 0.1, 0.4]
#       -        
# Nm:  1258342.0
# Nv:  315070.00000000006
# N total:  1573412.0
# Y1:  0.012299084829832252
# m safe:  True
# v safe:  True
# qs:  [0.0993612781077823, 0.06139346113812681, 0.10799261669330697, 0.04757324659463424]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  1665782.4148909424
# QBER:  0.05089783651567712



# ---> d =  8000
# x:  [0.6, 0.1, 0.4]
#       -        
# Nm:  805255.0
# Nv:  201379.99999999997
# N total:  1006635.0
# Y1:  0.007787270090497662
# m safe:  True
# v safe:  True
# qs:  [0.10110099915471271, 0.06261297047763685, 0.11012856385830411, 0.04693561093316429]
# ---- P_fail
# (0.0, 0.0, 0.0)
# ------ ES
# signal:  1061894.4148909424
# QBER:  0.051233455426985965

# ---> d =  9000
# x:  [0.6, 0.1, 0.4]
#       -        
# Nm:  513776.0
# Nv:  128886.0
# N total:  642662.0
# Y1:  0.004928187734842759
# m safe:  True
# v safe:  True
# qs:  [0.1028485470467531, 0.06414135615416962, 0.11230932613893929, 0.04608387249363588]
# ---- P_fail
# (0.0, 0.0, 0.0)
# Key rate:  -0.00012508012962152392
# ------ ES
# signal:  680647.4148909423
# QBER:  0.05170924393997868