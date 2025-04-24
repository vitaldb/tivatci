# tivatci Library Documentation

This library provides tools for pharmacokinetic (PK) and pharmacodynamic (PD) modeling, particularly for Target Controlled Infusion (TCI).

## Installation

```bash
pip install tivatci
```

## Usage

```python
from tivatci import Model, calc_lbm

# Example using the Schnider model for propofol
patient_age = 40
patient_sex = 'M'
patient_weight = 70  # kg
patient_height = 175 # cm

model = Model(name='schnider', age=patient_age, sex=patient_sex, weight=patient_weight, height=patient_height)
print(model)

# Simulate infusion to target 3 ug/mL effect-site concentration
target_concentration = 3 # ug/mL
simulation_duration = 600 # seconds
cts = [target_concentration] * simulation_duration # Target concentration array

results_df = model.run(cts, filename='schnider_simulation')
print(results_df.head())

# Calculate Lean Body Mass (LBM) using James model
lbm = calc_lbm(sex=patient_sex, weight=patient_weight, height=patient_height, model='james')
print(f"Lean Body Mass (James): {lbm:.2f} kg")
```

## API Reference

### `Model` Class

Represents a pharmacokinetic/pharmacodynamic model.

#### `__init__(self, name=None, age=30, sex='F', weight=60, height=160, v1=None, k10=0, k12=0, k13=0, k21=0, k31=0, v2=0, v3=0, q1=0, q2=0, q3=0, ke0=0)`

Initializes the model. Parameters can be set based on a predefined model `name` and patient demographics (`age`, `sex`, `weight`, `height`), or by providing specific PK/PD parameters directly (either rate constants `k` or clearances `q` and volumes `v`).

-   **`name`** (str, optional): Name of the predefined model (e.g., 'marsh', 'schnider', 'minto', 'eleveld', 'gepts', etc.).
-   **`age`** (int, optional): Patient age in years. Defaults to 30.
-   **`sex`** (str, optional): Patient sex ('M' or 'F'). Defaults to 'F'.
-   **`weight`** (float, optional): Patient weight in kg. Defaults to 60.
-   **`height`** (float, optional): Patient height in cm. Defaults to 160.
-   **`v1`, `k10`, `k12`, `k13`, `k21`, `k31`** (float, optional): Volume of central compartment and rate constants (min⁻¹).
-   **`v1`, `v2`, `v3`, `q1`, `q2`, `q3`** (float, optional): Volumes of distribution (L) and intercompartmental clearances (L/min).
-   **`ke0`** (float, optional): Effect-site elimination rate constant (min⁻¹).

#### `setq(self, v1=0, v2=0, v3=0, q1=0, q2=0, q3=0, ke0=0)`

Sets model parameters using volumes (`v1`, `v2`, `v3`) and clearances (`q1`, `q2`, `q3`). Converts these to rate constants (`k`) internally.

#### `setk(self, v1=0, k10=0, k12=0, k13=0, k21=0, k31=0, ke0=0)`

Sets model parameters directly using the volume of the central compartment (`v1`) and rate constants (`k10`, `k12`, `k13`, `k21`, `k31`).

#### `cp(self, tmax=None, dose=None, a=None)`

Calculates the plasma concentration over time.

-   **`tmax`** (int, optional): Maximum simulation time in seconds.
-   **`dose`** (list or float, optional): Infused amount per second. A single float is treated as a bolus at time 0.
-   **`a`** (np.array, optional): Initial amounts in compartments. Defaults to zeros.
-   **Returns**: `np.array` of plasma concentrations.

#### `ce(self, tmax=None, dose=None, a=None, ke0=None)`

Calculates the effect-site concentration over time.

-   **`tmax`** (int, optional): Maximum simulation time in seconds.
-   **`dose`** (list or float, optional): Infused amount per second.
-   **`a`** (np.array, optional): Initial amounts in compartments. Defaults to zeros.
-   **`ke0`** (float, optional): Effect-site elimination rate constant. Defaults to `self.ke0`.
-   **Returns**: `np.array` of effect-site concentrations.

#### `sim(self, tmax=None, dose=None, a=None, ke0=None)`

Simulates the amount of drug in each compartment over time using a discrete-time matrix approach (updates every second).

-   **`tmax`** (int, optional): Maximum simulation time in seconds. If `None`, simulates until peak effect-site concentration is reached after a bolus. Defaults to 9999 seconds if `dose` is not a single bolus.
-   **`dose`** (list or float, optional): Infused amount per second. Defaults to 0.
-   **`a`** (np.array, optional): Initial amounts in compartments `[a1, a2, a3, a4]`. Defaults to `[0, 0, 0, 0]`.
-   **`ke0`** (float, optional): Effect-site elimination rate constant. Defaults to `self.ke0`.
-   **Returns**: `np.array` of shape `(time_steps, 4)` representing drug amounts in compartments [central, peripheral1, peripheral2, effect-site].

#### `tpeak(self, ke0=None, prec=1)`

Estimates the time (in seconds) to reach the maximum effect-site concentration after a bolus dose.

-   **`ke0`** (float, optional): Effect-site elimination rate constant. Defaults to `self.ke0`.
-   **`prec`** (int, optional): Precision of the time calculation in seconds. Defaults to 1.
-   **Returns**: float, time to peak effect in seconds.

#### `recalc_ke0(self, tpeak)`

Finds the optimal `ke0` (min⁻¹) that matches a given time to peak effect. Uses `scipy.optimize.brentq`.

-   **`tpeak`** (float): The target time to peak effect in seconds.
-   **Returns**: float, the optimized `ke0` value.

#### `udf(self, plasma=False)`

Generates the Unit Disposition Function (UDF) for either the plasma or effect site after a 10-second unit infusion.

-   **`plasma`** (bool, optional): If `True`, returns plasma UDF. If `False`, returns effect-site UDF. Defaults to `False`.
-   **Returns**: `np.array` representing the concentration profile.

#### `tci(self, ct, a=None, plasma=False)`

Calculates the required infusion rate to reach a target concentration (`ct`) using the Shafer and Gregg algorithm, considering the current state (`a`).

-   **`ct`** (float): Target concentration (plasma or effect-site).
-   **`a`** (np.array, optional): Current amounts in compartments. Defaults to zeros.
-   **`plasma`** (bool, optional): If `True`, targets plasma concentration. If `False`, targets effect-site concentration. Defaults to `False`.
-   **Returns**: tuple `(rate, tpeak)`, where `rate` is the calculated infusion rate and `tpeak` is the predicted time until the target is reached or the infusion should be re-evaluated.

#### `run(self, cts, filename=None, maxrate=None, plasma=False)`

Runs a full TCI simulation for a given target concentration profile (`cts`).

-   **`cts`** (list or np.array): Array of target concentrations for each second of the simulation.
-   **`filename`** (str, optional): If provided, saves simulation results (Cp, Ce, Ct, Rate, Infused) to a CSV file and a plot of concentrations to a PNG file with this base name.
-   **`maxrate`** (float, optional): Maximum allowed infusion rate.
-   **`plasma`** (bool, optional): If `True`, targets plasma concentration. Defaults to `False` (targets effect-site).
-   **Returns**: `pd.DataFrame` containing columns: 'Ct', 'Cp', 'Ce', 'Rate', 'Infused'.

#### `__repr__(self)`

Returns a string representation of the model's parameters.

### Standalone Functions

#### `_sigmoid(x, e50, y)`

Calculates the sigmoid Emax model value. Used internally for maturation functions (e.g., in the Eleveld model).

-   **`x`**: Input value (e.g., age, weight).
-   **`e50`**: Value of `x` at which 50% of the maximum effect is achieved.
-   **`y`**: Hill coefficient (slope factor).
-   **Returns**: float, the sigmoid function output.

#### `calc_lbm(sex, weight, height, age=None, model='james')`

Calculates Lean Body Mass (LBM) based on patient demographics using various formulas.

-   **`sex`** (str): Patient sex ('M' or 'F').
-   **`weight`** (float): Patient weight in kg.
-   **`height`** (float): Patient height in cm.
-   **`age`** (int, optional): Patient age in years (required for some models like 'al-sallami').
-   **`model`** (str, optional): The LBM formula to use ('james', 'janmahasatian', 'devine', 'al-sallami'). Defaults to 'james'.
-   **Returns**: float, calculated Lean Body Mass in kg.
-   **Raises**: `ValueError` if an unsupported model is specified or if age is required but not provided.

---
*Documentation generated based on the source code.*
