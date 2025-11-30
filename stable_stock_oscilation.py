import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Parameters Dataset 2 from  paper 
p = {
    's1': 0.04, 's2': 0.05,
    'K1': 5,    'K2': 10,
    'm12': 0.01, 'm21': 0.02,
    'v1': 0.03, 'v2': 0.03,
    'a1': 0.1,  'a2': 0.3,
    'b1': 0.02, 
    'b2': 1.5,  
    'c1': 0.15,  
    'c2': 0.2,   
    'mu': 0.02
}

# Initial Conditions
y0 = [0.54, 3.3, 0.75] # X1, X2, Y

# Time Setup
t_span = [0, 20000]
t_eval = np.linspace(0, 20000, 20000)

#  System of DE
def prey_predator_system(t, y):
    X1, X2, Y = y
    
    # Functional Response
    r1 = X1 / (X1 + p['a2'] * X2 + p['b2'] * Y)
    r2 = X2 / (X2 + p['a1'] * X1 + p['b1'] * Y)
    
    dX1_dt = p['s1']*X1*(1 - X1/p['K1']) - p['m12']*X1*X2 - p['v1']*X1*Y*r1
    dX2_dt = p['s2']*X2*(1 - X2/p['K2']) - p['m21']*X1*X2 - p['v2']*X2*Y*r2
    dY_dt  = -p['mu']*Y + p['c1']*p['v1']*X1*Y*r1 + p['c2']*p['v2']*X2*Y*r2
    
    return [dX1_dt, dX2_dt, dY_dt]

# Solve
solution = solve_ivp(prey_predator_system, t_span, y0, t_eval=t_eval, method='RK45')

#   Data for Seaborn
df = pd.DataFrame(solution.y.T, columns=['Prey 1 (X1)', 'Prey 2 (X2)', 'Predator (Y)'])
df['Time'] = solution.t

 
df_long = df.melt('Time', var_name='Company', value_name='Relative Share Price')

# Plot
sns.set_theme(style="whitegrid", context="notebook", font_scale=1.1)
plt.figure(figsize=(12, 6))

# Plot  
sns.lineplot(
    data=df_long, 
    x='Time', 
    y='Relative Share Price', 
    hue='Company', 
    palette='deep',
    linewidth=2
)

 
plt.title(f"Stable Equilibrium (Dataset 2): $b_2={p['b2']}, c_1={p['c1']}, c_2={p['c2']}$", 
          fontsize=16, fontweight='bold', pad=20)
plt.ylabel("Population Sizes / Share Price")
plt.xlabel("Time, t")
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)

plt.tight_layout()
plt.show()