import tivatci

model = tivatci.Model('schnider', 80, 'M', 75, 172)
cts = [4] * 200 + [3] * 200 + [5] * 160 + [2] * 200 + [0] * 500
model.run(cts, 'result.csv')
