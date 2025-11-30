import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

 
# Fixed b2 = 3.0. 
p = {
    's1': 0.04, 's2': 0.05,
    'K1': 5,    'K2': 10,
    'm12': 0.01, 'm21': 0.02,
    'v1': 0.03, 'v2': 0.03,
    'a1': 0.1,  'a2': 0.3,
    'b1': 0.02, 'b2': 3.0,
    'c1': 0.1,
    'mu': 0.02
}

 
 
# We scan from 0.220 to 0.235 to see the fork open up
c2_values = np.linspace(0.220, 0.235, 300)

final_peaks = []
final_c2 = []

# Initial condition 
y0 = [0.54, 3.3, 0.75] 

# Sim settings
t_span = [0, 8000]   
t_eval = np.linspace(6000, 8000, 1000)  

 

# System of DE
def system(t, y, c2):
    X1, X2, Y = y
    
    denom1 = X1 + p['a2'] * X2 + p['b2'] * Y
    denom2 = X2 + p['a1'] * X1 + p['b1'] * Y
    
    # Functional Responses
    r1 = X1 / denom1 if denom1 > 1e-9 else 0
    r2 = X2 / denom2 if denom2 > 1e-9 else 0
    
    dX1 = p['s1']*X1*(1 - X1/p['K1']) - p['m12']*X1*X2 - p['v1']*X1*Y*r1
    dX2 = p['s2']*X2*(1 - X2/p['K2']) - p['m21']*X1*X2 - p['v2']*X2*Y*r2
    dY  = -p['mu']*Y + p['c1']*p['v1']*X1*Y*r1 + c2*p['v2']*X2*Y*r2
    
    return [dX1, dX2, dY]

 
for c2 in c2_values:
    sol = solve_ivp(
        lambda t, y: system(t, y, c2), 
        t_span, 
        y0, 
        t_eval=t_eval, 
        method='RK45', 
        rtol=1e-8, 
        atol=1e-10
    )
    
 
    x1_data = sol.y[0]
    
    
    if sol.y[2][-1] < 1e-4:
        continue

    # Peak Detection
    peaks, _ = find_peaks(x1_data)
    valleys, _ = find_peaks(-x1_data)
    
    if len(peaks) > 1:
        # OSCILLATION
 
        unique_peaks = np.unique(np.round(x1_data[peaks], 4))
        unique_valleys = np.unique(np.round(x1_data[valleys], 4))
        
        for val in unique_peaks:
            final_c2.append(c2)
            final_peaks.append(val)
        for val in unique_valleys:
            final_c2.append(c2)
            final_peaks.append(val)
    else:
        # STABLE 
        final_c2.append(c2)
        final_peaks.append(x1_data[-1])
    
 
    y0 = sol.y[:, -1]

# Plot
plt.figure(figsize=(12, 7))
plt.scatter(final_c2, final_peaks, s=2, c='black', alpha=0.6)

 

plt.axvline(x=0.228283, color='red', linestyle='--', linewidth=1.5, label='Theoretical HB Point (0.228)')

plt.title('Hopf Bifurcation in Stock Model (Dataset 2)', fontsize=14)
plt.xlabel('Parameter $c_2$ (Conversion Rate)', fontsize=12)
plt.ylabel('Prey 1 Share Price ($X_1$)', fontsize=12)
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

 