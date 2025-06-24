import matplotlib.pyplot as plt
import numpy as np

plt.close('all')

# ra = 20e6  t = 1
zs = [0, 10, 20, 30, 40, 50, 60]

up_sig_ES = [1936.6045719801918, 1831.6045719801918, 1785.6045719801918, 1289.6045719801918, 851.6045719801917, 474.60457198019174, 178.60457198019176]

up_qber_ES = [0.06699186443335407, 0.06749517822332543, 0.06773002206911477, 0.07104664820500173, 0.07626481075145483, 0.0861419345204327, 0.11342972351522757]

up_sig_decoy = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

up_qber_decoy = [0.25544301347831505, 0.2676749628929547, 0.2715571006355309, 0.3757066207043998, 0.7312710632178389, -27.161629652379215, -0.5535480186139646]

down_sig_ES = [19450.60457198019, 17791.60457198019, 17399.60457198019, 14899.604571980191, 11253.604571980191, 6624.604571980191, 3260.604571980192]

down_qber_ES = [0.055204456128256356, 0.05544485483573655, 0.05550667004468906, 0.05595724418286727, 0.056869941990569856, 0.05900128158472092, 0.06295917583325711]

down_sig_decoy = [5303.0, 4868.000000000001, 4746.0, 4073.9999999999995, 3103.9999999999995, 1831.0000000000002, 900.0]

down_qber_decoy = [0.09115955133598962, 0.09296771972324504, 0.09378982764325093, 0.0976190683267772, 0.10532073001475753, 0.12648068939546578, 0.17857880801991163]




# ra = 100e6  t = 1
up_sig_ES2 = [9781.022859900959, 9247.022859900959, 8880.022859900959, 6412.022859900959, 4227.022859900959, 2394.022859900959, 891.0228599009588]

up_qber_ES2 = [0.057890263416312454, 0.05813486007794624, 0.05831636240282545, 0.059945524476186174, 0.0625586473979712, 0.06743751709719947, 0.08204055916011777]

up_sig_decoy2 = [2682.0, 2528.0, 2474.0, 1764.0000000000002, 1171.0, 0.0, 0.0]

up_qber_decoy2 = [0.11273108307926478, 0.11523938702575068, 0.11566383767307457, 0.13297890356509057, 0.1643726572863023, 0.24142407828752616, 1.0486628290565478]

down_sig_ES2 = [96391.02285990096, 89220.02285990096, 86947.02285990096, 74361.02285990096, 56400.02285990096, 33260.02285990096, 16073.022859900959]

down_qber_ES2 = [0.052371627817303, 0.052467696544344994, 0.05250064777315071, 0.05271021397900766, 0.05312611032421931, 0.05411440526192896, 0.05603981517262192]

down_sig_decoy2 = [49005.0, 44968.0, 44242.0, 37837.0, 29141.000000000007, 8957.0, 4442.0]

down_qber_decoy2 = [0.07074966359653555, 0.07134274683277632, 0.07163015156145072, 0.07277861210374006, 0.07559506932759211, 0.08273068420982377, 0.09650143621553277]



zs = 90 - np.array(zs)

fig = plt.figure()
fig.set_size_inches(18.5*.3, 10.5*.3)

plt.plot(zs, up_sig_ES, linewidth=2, label='Uplink - Entangle source')
# plt.plot(zs, up_sig_decoy, linewidth=2, label='Uplink - Entangle source')
plt.plot(zs, down_sig_ES, linewidth=2, label='Downlink - Decoy')
plt.plot(zs, down_sig_decoy, '--', linewidth=2, label='Downlink - Decoy')

plt.xlabel('Elevation (deg)',fontsize=13)
plt.ylabel('Signal', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)
plt.legend()

fig = plt.figure()
fig.set_size_inches(18.5*.3, 10.5*.3)

plt.plot(zs, up_qber_ES, linewidth=2, label='Uplink - Entangled source')
# plt.plot(zs, up_qber_decoy, linewidth=2, label='Uplink - Entangle source')
plt.plot(zs, down_qber_ES, linewidth=2, label='Downlink - Entangled source')
plt.plot(zs, down_qber_decoy, '--', linewidth=2, label='Downlink - Decoy')

plt.xlabel('Elevation (deg)',fontsize=13)
plt.ylabel('QBER', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)
plt.legend()




fig = plt.figure()
fig.set_size_inches(18.5*.3, 10.5*.3)

plt.plot(zs, up_sig_ES2, linewidth=2, label='Uplink - Entangle source')
plt.plot(zs[:-2], up_sig_decoy2[:-2], '--', linewidth=2, label='Uplink - Decoy')
plt.plot(zs, down_sig_ES2, linewidth=2, label='Downlink - Entangled source')
plt.plot(zs, down_sig_decoy2, '--', linewidth=2, label='Downlink - Decoy')

plt.xlabel('Elevation (deg)',fontsize=13)
plt.ylabel('Signal', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)
plt.legend()

fig = plt.figure()
fig.set_size_inches(18.5*.3, 10.5*.3)

plt.plot(zs, up_qber_ES2, linewidth=2, label='Uplink - Entangle source')
plt.plot(zs[:-2], up_qber_decoy2[:-2], '--', linewidth=2, label='Uplink - Decoy')
plt.plot(zs, down_qber_ES2, linewidth=2, label='Downlink - Entangled source')
plt.plot(zs, down_qber_decoy2, '--', linewidth=2, label='Downlink - Decoy')

plt.xlabel('Elevation (deg)',fontsize=13)
plt.ylabel('QBER', fontsize=13)

plt.yticks(fontsize=13)
plt.xticks(fontsize=13)
plt.legend()


