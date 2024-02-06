import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class Model:
    def __init__(self, name, age, sex, wt, ht):
        self.set_model(name, age, sex, wt, ht)

        # convert units from min to sec
        self.k10 /= 60.0
        self.k12 /= 60.0
        self.k21 /= 60.0
        self.k13 /= 60.0
        self.k31 /= 60.0
        self.ke0 /= 60.0
        self.v4 = self.v1 / 1000
        self.k14 = self.ke0 * self.v4 / self.v1

        # generate update matrix
        self.k = np.array([
            [1 - self.k10 - self.k12 - self.k13, self.k21, self.k31, 0],
            [self.k12, 1 - self.k21, 0, 0],
            [self.k13, 0, 1 - self.k31, 0],
            [self.k14, 0, 0, 1 - self.ke0]
        ])

        # generate udf using the numerical estimation
        self.reset()
        self.e_udf = np.concatenate((self.ce(10, 1, True), self.ce(2700 - 10)))
        self.reset()        
    
    def set_model(self, name, age, sex, wt, ht):
        self.v1 = 0
        self.k10 = 0
        self.k12 = 0
        self.k13 = 0
        self.k21 = 0
        self.k31 = 0
        self.ke0 = 0

        name = name.lower()
        if name == 'marsh':
            self.v1 = 0.228 * wt
            self.k10 = 0.119
            self.k12 = 0.114
            self.k13 = 0.0419
            self.k21 = 0.055
            self.k31 = 0.0033
            self.ke0 = 0.26  # diprifusor
        elif name == "modified marsh":
            self.v1 = 0.228 * wt
            self.k10 = 0.119
            self.k12 = 0.114
            self.k13 = 0.0419
            self.k21 = 0.055
            self.k31 = 0.0033
            self.ke0 = 1.2195  # stanpump, orchestra
        elif name == "schnider":
            lbm = james(sex, wt, ht)
            self.v1 = 4.27
            v2 = 18.9 - 0.391 * (age - 53)
            v3 = 238
            cl1 = 1.89 + 0.0456 * (wt - 77) - 0.0681 * (lbm - 59) + 0.0264 * (ht - 177)
            cl2 = 1.29 - 0.024 * (age - 53)
            cl3 = 0.836
            self.k10 = cl1 / self.v1
            self.k12 = cl2 / self.v1
            self.k13 = cl3 / self.v1
            self.k21 = cl2 / v2
            self.k31 = cl3 / v3
            self.ke0 = 0.456
        elif name == "paedfusor":
            if 1 <= age < 13:
                self.v1 = 0.4584 * wt
                self.k10 = 0.1527 * wt ** -0.3
            elif age <= 13:
                self.v1 = 0.4 * wt
                self.k10 = 0.0678
            elif age <= 14:
                self.v1 = 0.342 * wt
                self.k10 = 0.0792
            elif age <= 15:
                self.v1 = 0.284 * wt
                self.k10 = 0.0954
            elif age <= 16:
                self.v1 = 0.22857 * wt
                self.k10 = 0.119
            else:
                self.v1 = None
            self.k12 = 0.114
            self.k13 = 0.0419
            self.k21 = 0.055
            self.k31 = 0.0033
            self.ke0 = 0.26  # from diprifusor (for adults)
            self.ke0 = 0.91  # Munoz et al Anesthesiology 2004:101(6)
        elif name == "kataria":  # Kataria et al. Anesthesiology 199480:104
            self.v1 = 0.41 * wt
            v2 = 0.78 * wt + 3.1 * age - 15.5
            v3 = 6.9 * wt
            cl1 = 0.035 * wt
            cl2 = 0.077 * wt
            cl3 = 0.026 * wt
            self.k10 = cl1 / self.v1
            self.k12 = cl2 / self.v1
            self.k13 = cl3 / self.v1
            self.k21 = cl2 / v2
            self.k31 = cl3 / v3
            self.ke0 = 0.41  # Munoz et al Anesthesiology 2004:101(6)
        elif name == "kim":
            self.v1 = 1.69
            v2 = 27.2 + 0.93 * (wt - 25)
            cl1 = 0.89 * (wt / 23.6) ** 0.97
            cl2 = 1.3
            self.k10 = cl1 / self.v1
            self.k12 = cl2 / self.v1
            self.k13 = 0
            self.k21 = cl2 / v2
            self.k31 = 0
        elif name == "minto":
            lbm = james(sex, wt, ht)
            self.v1 = 5.1 - 0.0201 * (age - 40) + 0.072 * (lbm - 55)
            v2 = 9.82 - 0.0811 * (age - 40) + 0.108 * (lbm - 55)
            v3 = 5.42
            cl1 = 2.6 - 0.0162 * (age - 40) + 0.0191 * (lbm - 55)
            cl2 = 2.05 - 0.0301 * (age - 40)
            cl3 = 0.076 - 0.00113 * (age - 40)
            self.k10 = cl1 / self.v1
            self.k12 = cl2 / self.v1
            self.k13 = cl3 / self.v1
            self.k21 = cl2 / v2
            self.k31 = cl3 / v3
            self.ke0 = 0.595 - 0.007 * (age - 40)

    def reset(self):
        '''drug amount in the compartments'''
        self.a = np.zeros(4)
    
    def sim(self, t, rate=0, update=False):
        '''simulate the movement of drug amount in the compartments
        t: time in seconds or an array of time
        rate: infusion amount/sec
        update: update the internal state if True
        returns: the estimated amount of drugs in the compartments in the shape (t, 4)
        '''
        ret = []
        a = self.a
        for _ in range(t):
            a = np.matmul(self.k, a)  # natural decay
            a[0] += rate
            ret.append(a)
        if update:
            self.a = a
        return np.array(ret)
    
    def ce(self, t, rate=0, update=False):
        return self.sim(t, rate, update)[:,3] / self.v4

    def tci(self, ct):  # shafer and greg algorithm
        '''calculate the infusion rate to achieve the desired effect site concentration'''
        tpeak = np.argmax(self.e_udf)  # initial tpeak
        while True:
            bes = self.ce(tpeak+1)  # natural decay
            if bes[0] > ct: # if the current effect site concentration is higher than the target
                if bes[-1] < ct: # wait until the effect site concentration is lower than the target
                    return 0, np.where(bes < ct)[0][0]
                return 0, tpeak
            rate = (ct - bes[-1]) / self.e_udf[tpeak]
            ces = bes + self.e_udf[:len(bes)] * rate
            new_tpeak = np.argmax(ces)
            if new_tpeak == tpeak or abs(ces[new_tpeak] - ct) <= ct * 0.001:
                break
            tpeak = new_tpeak
        return rate, tpeak
    
    def run(self, cts, filename=None):
        last_ct = 0
        wait_until = 0
        infuse_until = 0
        ces = []
        cps = []
        rates = []
        for i in range(len(cts)):
            ct = cts[i]
            if i >= wait_until or ct != last_ct:
                last_ct = ct
                rate, tpeak = self.tci(cts[i])
                infuse_until = i + 10
                #print(f'{i}sec @ rate: {rate}, tpeak: {tpeak}')
                wait_until = i + tpeak
            if i >= infuse_until:
                rate = 0
            a = self.sim(1, rate, True)[-1]
            rates.append(rate)
            cps.append(a[0] / self.v1)
            ces.append(a[3] / self.v4)

        plt.figure(figsize=(20, 5))
        plt.plot(cps, color='red', label='Cp')
        plt.plot(cts, color='blue', label='Ct')
        plt.plot(ces, color='green', label='Ce')
        plt.legend()
        if filename:
            plt.savefig(filename + '.png')
            pd.DataFrame({'Cp': cps, 'Ct': cts, 'Ce': ces, 'Rate': rates}).to_csv(filename, index=False)    
        else:
            plt.show()


def james(sex, wt, ht):
    if sex == "M":
        return 1.1 * wt - 128 * (wt / ht) ** 2
    elif sex == "F":
        return 1.07 * wt - 148 * (wt / ht) ** 2
    raise ValueError


if __name__ == '__main__':
    model = Model('schnider', 80, 'M', 75, 172)
    cts = np.concatenate([
        np.full(200, 4),
        np.full(200, 3),
        np.full(160, 5),
        np.full(200, 2),
        np.full(300, 0),
        ])
    model.run(cts, 'result.csv')
