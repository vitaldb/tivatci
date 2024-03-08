import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize as optimize

class Model:
    def __init__(self, name=None, age=30, sex='F', weight=60, height=160, v1=None, k10=0, k12=0, k13=0, k21=0, k31=0, v2=0, v3=0, q1=0, q2=0, q3=0, ke0=0):
        '''
        set the model parameters based on the name of the model and the patient
        time constants are always in min^-1
        '''
        if name is not None:
            name = name.lower()
        # either k or name should be provided
        if v1 and ke0:
            if k10 or k12 or k13 or k21 or k31:
                self.setk(v1, k10, k12, k13, k21, k31, ke0)
            elif v2 or v3 or q1 or q2 or q3:
                self.setq(v1, v2, v3, q1, q2, q3, ke0)
        elif name == 'marsh':  
            # Marsh et al. Pharmacokinetic model driven infusion of propofol in children. Br J Anaesth 1991;67:41-8
            self.setk(0.228 * weight, 0.119, 0.114, 0.0419, 0.055, 0.0033, 0.26) # diprifusor
        elif name == 'modified marsh':
            self.setk(0.228 * weight, 0.119, 0.114, 0.0419, 0.055, 0.0033, 1.2195)  # stanpump, orchestra
        elif name == 'schnider':
            lbm = calc_lbm(sex, weight, height)
            self.setq(4.27, 18.9 - 0.391 * (age - 53), 238,
                1.89 + 0.0456 * (weight - 77) - 0.0681 * (lbm - 59) + 0.0264 * (height - 177),
                1.29 - 0.024 * (age - 53), 0.836, 
                0.456)
        elif name == 'paedfusor':
            if 1 <= age < 13:
                v1 = 0.4584 * weight
                k10 = 0.1527 * weight ** -0.3
            elif age <= 13:
                v1 = 0.4 * weight
                k10 = 0.0678
            elif age <= 14:
                v1 = 0.342 * weight
                k10 = 0.0792
            elif age <= 15:
                v1 = 0.284 * weight
                k10 = 0.0954
            elif age <= 16:
                v1 = 0.22857 * weight
                k10 = 0.119
            else:
                raise ValueError
            ke0 = 0.26  # from diprifusor (for adults)
            ke0 = 0.91  # Munoz et al Anesthesiology 2004:101(6)
            self.setk(v1, k10, 0.114, 0.0419, 0.055, 0.0033, ke0) 
        elif name == 'kataria':  
            # Kataria et al. The pharmacokinetics of propofol in children using three different data analysis approaches. Anesthesiology 1994;80:104
            self.setq(0.41 * weight, 0.78 * weight + 3.1 * age - 15.5, 6.9 * weight,
                0.035 * weight, 0.077 * weight, 0.026 * weight,
                0.41) # Munoz et al Anesthesiology 2004:101(6)
        elif name == 'choi':  # propofol
            # Choi et al. Population pharmacokinetic and pharmacodynamic model of propofol externally validated in children. J Pharmacokinet Pharmacodyn. 2015 Apr;42(2):163-77
            self.setq(1.69, 27.2 + 0.93 * (weight - 25), 0,
                0.89 * (weight / 23.6) ** 0.97, 1.3, 0, 
                0.371)
        elif name == 'eleveld':
            t1 = 6.28 # v1 ref
            t2 = 25.5 # v2 ref
            t3 = 273 # v3 ref
            t4 = 1.79 # cl ref
            t5 = 1.75 # q2 ref
            t6 = 1.11 # q3 ref
            t8 = 42.3 # cl maturation e50
            t9 = 9.06 # cl maturation slope
            t10 = -0.0156 # smaller v2 with age
            t11 = -0.00286 # lower cl with age
            t12 = 33.6 # v1 sigmoid e50=33.6 kg
            t13 = -0.0138 # smaller v3 with age
            t14 = 68.3 # maturation of q3
            t15 = 2.1 # cl ref (female)
            t16 = 1.3 # higher q2 for maturation of q3

            # reference
            age_ref = 35
            weight_ref = 70
            ht_ref = 170
            with_opioids = True

            # age(yrs) to pma(wks)
            pma = age * 52 + 40
            pma_ref = age_ref * 52 + 40

            cl3_mat = _sigmoid(pma, t14, 1)
            cl3_mat_ref = _sigmoid(pma_ref, t14, 1)

            # lean body mass
            ffm_ref = calc_lbm('M', weight_ref, ht_ref, age_ref, 'al-sallami')
            ffm = calc_lbm(sex, weight, height, age, 'al-sallami')

            v1 = t1 * (_sigmoid(weight, t12, 1) / _sigmoid(weight_ref, t12, 1))
            v2 = t2 * (weight / weight_ref) * np.exp(t10 * (age - age_ref))
            v3 = t3 * ffm / ffm_ref
            if with_opioids:
                v3 *= np.exp(t13 * age)
            cl1 = (t4 if sex == 'M' else t15) * ((weight / weight_ref) ** 0.75)
            cl1 *= _sigmoid(pma, t8, t9) / _sigmoid(pma_ref, t8, t9)  # age maturation
            if with_opioids:
                cl1 *= np.exp(t11 * age)
            cl2 = t5 * (v2 / t2) ** 0.75 * (1 + t16 * (1 - cl3_mat))
            cl3 = t6 * (v3 / t3) ** 0.75 * (cl3_mat / cl3_mat_ref)
            ke0 = 0.146 * ((weight / weight_ref) ** -0.25)

            self.setq(v1, v2, v3, cl1, cl2, cl3, ke0)
        elif name == 'minto':
            # Minto et al. Influence of age and gender on the pharmacokinetics and pharmacodynamicsof remifentanil. I. Model development. Anesthesiology, 86:10–23, 1997.
            lbm = calc_lbm(sex, weight, height, 'james')
            self.setq(5.1 - 0.0201 * (age - 40) + 0.072 * (lbm - 55),
                9.82 - 0.0811 * (age - 40) + 0.108 * (lbm - 55),
                5.42,
                2.6 - 0.0162 * (age - 40) + 0.0191 * (lbm - 55),
                2.05 - 0.0301 * (age - 40),
                0.076 - 0.00113 * (age - 40),
                0.595 - 0.007 * (age - 40))
        elif name == 'kim':  # pk model of remifentanil
            lbm = calc_lbm(sex, weight, height, 'janmahasatian')
            self.setq(4.76 * (weight / 74.5) ** 0.658,
                8.4 * (lbm / 52.3) ** 0.573 - 0.0936 * (age-37),
                4 - 0.0477 * (age-37),
                2.77 * (weight / 74.5) ** 0.336 - 0.0149 * (age-37),
                1.94 - 0.0280 * (age - 37),
                0.197,
                0.595 - 0.007 * (age - 40))  # minto's ke0
        elif name == 'schuttler': # remimazolam
            # Schüttler et al. Pharmacokinetics and Pharmacodynamics of Remimazolam (CNS 7056) after Continuous Infusion in Healthy Male Volunteers: Part I. Pharmacokinetics and Clinical Pharmacodynamics. Anesthesiology. 2020 Apr;132(4):636-651.
            self.setq(4.7 * weight / 75, 14.5, 15.5, 
                1.14, 1.04, 0.19, 
                0.27)
        elif name == 'schmith':
            # Zhou et al. Population pharmacokinetic/pharmacodynamic modeling for remimazolam in the induction and maintenance of general anesthesia in healthy subjects and in surgical subjects. J Clin Anesth. 2020
            v1 = 2.92 * weight / 70
            v2 = 19.1 * weight / 70
            # if asa3:
            #     v1 *= 1 - 0.56
            #     v2 *= 1.22
            v3 = 9.81 * weight / 70
            cl1 = 61.6 * weight / 70 * (1.11 if sex == 'F' else 1) / 60
            cl2 = 22.9 * weight / 70 / 60
            cl3 = 69.6 * weight / 70 / 60
            ke0 = 8.08 / 60
            bmi = weight / (height / 100) ** 2
            if bmi > 25:
                ke0 *= 1.17
            # if asian: self.ke0 *= 1 - 0.48            
            self.setq(v1, v2, v3, cl1, cl2, cl3, ke0)
        elif name == 'wierda':  # rocuronium
            self.setk(45 * weight / 1000, 0.1, 0.21, 0.028, 0.13, 0.01, 0.168)
        elif name == 'cooper':  # rocuronium
            self.setk(38.5 * weight / 1000, 0.119, 0.259, 0.06, 0.163, 0.012, 0.168)
        elif name == 'saldien':  # rocuronium
            self.setk(35.6 * weight / 1000, 0.126, 0.209, 0.05, 0.163, 0.015, 0.168)
        elif name in ('dehaes', 'de haes'):  # rocuronium
            self.setk(42.0 * weight / 1000, 0.0762, 0.124, 0.0214, 0.13, 0.013, 0.15)
        elif name == 'gepts': # sufentanil
            self.setq(14.3, 63.1, 261.6, 0.92, 1.55, 0.33, ke0) # 10.1097/00000542-199512000-00010
            
            #self.ke0 = self.recalc_ke0(5.6 * 60) # 0.176/min = 10.1097/00000542-199101000-00010 (Shafer/Varvel 1991)
            # the ke0 value above is wrong based on the https://blog.tivatrainer.com/sufentanil-another-mistake-in-the-implementation-in-commercial-target-controlled-infusion-systems

            #self.ke0 = 0.119  # this value is from the stanpump code, but it is not clear where it comes from
            self.ke0 = 0.112  # Implementation of Gepts model in Fresenius Base Primea, Scott, Cooke, & Stanski, 1991, K for half life = 6.2 min
        elif name == 'hannivoort': # dexemedetomidine
            # 10.1097/ALN.0000000000000740
            self.setq(1.78 * (weight / 70), 30.3 * (weight / 70), 52.0 * (weight / 70),
                      0.686 * (weight / 70) ** 0.75, 2.98 * (weight / 70) ** 0.75, 0.602 * (weight / 70) ** 0.75)
            self.ke0 = 0.120 # 10.1093/bja/aex085 (Colin et al. 2017)
        elif name == 'scott': # fentanyl 
            # # k = 0.693 / hl
            # # vd = cl * hl / 0.693
            # hl2 = 1.0 # min
            # hl3 = 18.5 # min
            # q1 = 0.574 # L/min
            # q2 = 4.005 # L/min
            # q3 = 1.952 # L/min
            # v1 = 12.7
            # v2 = q2 * hl2 / 0.693 # = np.log(2)
            # v3 = q3 * hl3 / 0.693 # = np.log(2)
            # self.setq(v1, v2, v3, q1, q2, q3) # 1987, J Pharmacol Exp Ther 1987 Jan;240(1):159-166
            
            self.setk(12.7, .056, .373, .180, .096, .0077) # from stanpump code
            # = self.setq(12.7, 49.34, 296.88, 0.71, 4.74, 2.29) # from pkpdtools.com
            self.ke0 = 0.693 / 4.7 # 4.7 min half life
        # elif name == 'shafer': # fentanyl (unweighted)
        #     self.setq(6.09, 28.1, 228, .504, 2.87, 1.37) # Anesthesiology. 1990;73:1091-102. 10.1097/00000542-199012000-00005
        #     self.ke0 = self.recalc_ke0(3.6 * 60) # 10.1097/00000542-199101000-00010 (Shafer/Varvel 1991)
        # elif name in ('mcclain', 'hug'): # fentanyl
        #     self.setk(0.356 * weight, 0.041, 0.185, 0.141, 0.103, 0.020)
        #     self.ke0 = self.recalc_ke0(3.6 * 60) # 10.1097/00000542-199101000-00010 (Shafer/Varvel 1991)
        else:
            raise ValueError('unsupported model')

    def setq(self, v1=0, v2=0, v3=0, q1=0, q2=0, q3=0, ke0=0):
        if v2 == 0:
            self.setk(v1, q1 / v1, q2 / v1, q3 / v1, 0, 0, ke0)
        elif v3 == 0:
            self.setk(v1, q1 / v1, q2 / v1, q3 / v1, q2 / v2, 0, ke0)
        else:
            self.setk(v1, q1 / v1, q2 / v1, q3 / v1, q2 / v2, q3 / v3, ke0)

    v1_v4 = 1000
    def setk(self, v1=0, k10=0, k12=0, k13=0, k21=0, k31=0, ke0=0):
        '''
        set PK/PD parameters of the model
        '''
        self.v1 = v1
        self.k10 = k10
        self.k12 = k12
        self.k13 = k13
        self.k21 = k21
        self.k31 = k31
        self.ke0 = ke0

    def cp(self, tmax=None, dose=None, a=None):
        '''
        calculate the plasma concentration
        '''
        return self.sim(tmax, dose, a)[:, 0] / self.v1

    def ce(self, tmax=None, dose=None, a=None, ke0=None):
        '''
        calculate the effect site concentration
        '''
        return self.sim(tmax, dose, a, ke0)[:, 3] / self.v1 * self.v1_v4
    
    def sim(self, tmax=None, dose=None, a=None, ke0=None):
        '''
        simulate the movement of drug amount in the compartments
        tmax: maximum time in seconds for simulation (if None, do while the maximum Ce is reached)
        dose: infused amount at every second ([1] for bolus at time 0, [1] * 10 for infusion for 10 sec, etc.)
        ke0: the elimination rate from the effect site (if None, use the self.ke0)
        a: initial state of the compartments
        returns: the estimated amount of drugs in the compartments in the shape (tmax, 4)
        '''
        ret = []

        if dose is None:
            dose = 0
        if np.isscalar(dose):
            dose = [dose]

        if a is None:
            a = np.zeros(4)

        last_a4 = 0
        if tmax is None:
            tmax = 9999

        if ke0 is None:
            ke0 = self.ke0

        # generate update matrix
        k = np.array([
            [1 - (self.k10 + self.k12 + self.k13) / 60, self.k21 / 60, self.k31 / 60, 0],
            [self.k12 / 60, 1 - self.k21 / 60, 0, 0],
            [self.k13 / 60, 0, 1 - self.k31 / 60, 0],
            [ke0 / self.v1_v4 / 60, 0, 0, 1 - ke0 / 60]
        ])

        # simulation
        for i in range(tmax):
            a = np.matmul(k, a)  # natural decay
            if i < len(dose):
                a[0] += dose[i]  # infusion
            if tmax == 9999:  # when user wants do until maximum Ce is reached
                if a[3] < last_a4:
                    break
                last_a4 = a[3]
            ret.append(a)

        return np.array(ret)
    
    def tpeak(self, ke0=None, prec=1):
        '''
        estimate the time to reach the maximum effect site concentration
        ke0: the elimination rate from the effect site (if None, use the self.ke0)
        prec: the precision of the time in seconds
        '''
        if ke0 is None:
            ke0 = self.ke0
        model = Model(v1=self.v1, k10=self.k10 * prec, k12=self.k12 * prec, k13=self.k13 * prec, k21=self.k21 * prec, k31=self.k31 * prec, ke0=ke0 * prec)
        tpeak = len(model.ce(dose=1, ke0=ke0 * prec)) * prec
        return tpeak
    
    def recalc_ke0(self, tpeak):
        '''
        find optimal ke0 that minimize the difference between the estimated tpeak and the actual tpeak
        tpeak: the actual time to reach the maximum effect in seconds
        returns: the optimal ke0 in /min
        '''
        return optimize.brentq(lambda ke0, tpeak: self.tpeak(ke0, prec=0.1) - tpeak, a=1e-5, b=100, args=(tpeak))
    
    def udf(self, plasma=False):
        '''
        generate the unit disposition function (UDF) for the plasma or the effect site
        '''
        if plasma:
            return self.cp(10, dose=[1]*10) # always maximum at 10 sec
        else:  # effect site
            return self.ce(dose=[1]*10)

    def tci(self, ct, a=None, plasma=False):
        '''
        calculate the infusion rate to achieve the desired effect site concentration using the Shafer and Greg algorithm
        a: the initial state of the compartments
        ct: target concentration
        '''
        if a is None:
            a = np.zeros(4)
        
        # generate udf
        udf = self.udf(plasma=plasma)

        tpeak = len(udf)-1  # initial tpeak
        while True:
            # natural decay
            if plasma:
                b = self.cp(tpeak + 1, a=a)
            else:  # effect site
                b = self.ce(tpeak + 1, a=a)

            if b[0] > ct: # if the current effect site concentration is higher than the target
                if b[-1] < ct: # wait until the effect site concentration is lower than the target
                    return 0, np.where(b < ct)[0][0]
                return 0, tpeak
            
            rate = (ct - b[-1]) / udf[tpeak]
            c = b + udf[:len(b)] * rate
            new_tpeak = np.argmax(c)
            if new_tpeak == tpeak or abs(c[new_tpeak] - ct) <= ct * 0.001:
                break
            tpeak = new_tpeak
        return rate, tpeak
    
    def run(self, cts, filename=None, maxrate=None, plasma=False):
        '''
        simulate the movement of drug amount in the compartments
        returns: DataFrame of Ct, Cp, Ce, and Rate
        '''
        last_ct = 0
        wait_until = 0
        infuse_until = 0
        ces = []
        cps = []
        rates = []
        a = np.zeros(4) # initial state
        for i in range(len(cts)):
            ct = cts[i]
            if i >= wait_until or ct != last_ct:
                last_ct = ct
                rate, tpeak = self.tci(cts[i], a, plasma=plasma)
                infuse_until = i + 10
                if maxrate and rate > maxrate:
                    rate = maxrate
                    wait_until = i + 10
                else:
                    wait_until = i + tpeak + 1
            if i >= infuse_until:
                rate = 0
            a = self.sim(1, dose=rate, a=a)[-1]
            rates.append(rate)
            cps.append(a[0] / self.v1)
            ces.append(a[3] / self.v1 * self.v1_v4)
        
        df = pd.DataFrame({'Ct': cts, 'Cp': cps, 'Ce': ces, 'Rate': rates, 'Infused': np.cumsum(rates)})
        if filename:
            plt.figure(figsize=(20, 5))
            plt.plot(cps, color='red', label='Cp')
            plt.plot(cts, color='blue', label='Ct')
            plt.plot(ces, color='green', label='Ce')
            plt.legend()
            plt.savefig(filename + '.png')
            df.to_csv(filename, index=False)
        
        return df
    
    def __repr__(self):
        return f'Model(v1={self.v1:.2f}, k10={self.k10:.4f}, k12={self.k12:.4f}, k13={self.k13:.4f}, k21={self.k21:.4f}, k31={self.k31:.4f}, ke0={self.ke0:.4f})'

def _sigmoid(x, e50, y):
    return (x ** y) / ((x ** y) + (e50 ** y))


def calc_lbm(sex, weight, height, age=None, model='james'):
    if model == 'james':
        if sex == 'M':
            return 1.1 * weight - 128 * (weight / height) ** 2
        elif sex == 'F':
            return 1.07 * weight - 148 * (weight / height) ** 2
    elif model == 'janmahasatian':
        # Janmahasatian et al. Quantification of lean bodyweight. Clin Pharmacokinet. 2005;44(10):1051-65.
        bmi = weight / (height / 100) ** 2
        if sex == 'M':
            return 9270 * weight / (6680 + 216 * bmi)
        if sex == 'F':
            return 9270 * weight / (8780 + 244 * bmi)
    elif model == 'devine':
        # Devine BJ. Gentamicin therapy. Drug Intell Clin Pharm. 1974;8:650–655.
        if sex == 'M':
            return 50 + 0.91 * (height - 152.4)
        elif sex == 'F':
            return 45.5 + 0.91 * (height - 152.4)
    elif model == 'al-sallami':
        if age is None:
            raise ValueError
        bmi = weight / (height / 100) ** 2
        if sex == 'M':
            return (0.88 + ((1 - 0.88) / (1 + (age / 13.4) ** -12.7))) * 42.92 * weight / (30.93 + bmi)
        elif sex == 'F':
            return (1.11 + ((1 - 1.11) / (1 + (age / 7.1) ** -1.1))) * 37.99 * weight / (35.98 + bmi)
    raise ValueError

if __name__ == '__main__':
    print(Model('gepts'))
    # plt.figure(figsize=(20, 5))
    # plt.plot(model.udf(True), label='plasma')
    # plt.show()
