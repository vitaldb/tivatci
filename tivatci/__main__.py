import tivatci
import matplotlib.pyplot as plt

# plt.plot([lbm('M', x, 170, 'james') for x in range(30, 200)], label='james')
# plt.plot([lbm('M', x, 170, 'janmahasatian') for x in range(30, 200)], label='janmahasatian')
# plt.plot([lbm('M', x, 170, 'devine') for x in range(30, 200)], label='devine')
# plt.legend()
# plt.show()
# exit(0)

# model = Model('schnider', 80, 'M', 75, 172)
# print(f'tpeak at {model.tpeak()} sec') # original tpeak

# model.ke0 = model.recalculate_ke0(100)
# print(model.ke0)
# print(model.tpeak())
# exit(0)

# Compare models
dfs = []
names = ['modified marsh', 'schnider', 'eleveld']
cts = [4] * 200 + [3] * 200 + [5] * 160 + [2] * 200 + [0] * 500
names = ['schuttler', 'schmith']
cts = [1] * 200 + [1.5] * 200 + [1.2] * 160 + [1] * 200 + [0] * 500
conc_unit = 'ug/ml'
dose_unit = 'mg'
styles = ['solid', 'dashed', 'dotted', 'dashdot']
for name in names:
    model = tivatci.Model(name, 80, 'M', 75, 172)
    dfs.append(model.run(cts))

fig, ax1 = plt.subplots(figsize=(20, 5))
ax1.set_xlabel('time (s)')
ax1.plot(dfs[0]['Ct'], color='blue', label='Ct')
for i in range(len(names)):
    ax1.plot(dfs[i]['Cp'], linestyle=styles[i], color='red', label=f'Cp ({names[i]})')
    ax1.plot(dfs[i]['Ce'], linestyle=styles[i], color='green', label=f'Ce ({names[i]})')
ax1.set_ylabel(f'Concentration ({conc_unit})')
ax1.legend(loc='upper left')
ax2 = ax1.twinx()
for i in range(len(names)):
    ax2.plot(dfs[i]['Infused'], linestyle=styles[i], color='gray', label=f'Infused ({names[i]})')
ax2.set_ylabel(f'Infused ({dose_unit})')
ax2.legend(loc='upper right')
plt.show()
