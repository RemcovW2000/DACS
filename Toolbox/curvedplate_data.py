import numpy as np
from scipy import interpolate
import plotly.graph_objects as go

# Example data points (replace these with your actual digitized data)
data = {
    500: {'Zb': np.array([1, 18.33, 1e4]), 'k': np.array([3.949, 4.1118, 2.532e3])},
    750: {'Zb': np.array([1, 19.95, 1e4]), 'k': np.array([4.1753, 4.52, 2.359e3])},
    1000: {'Zb': np.array([1, 24.92, 1e4]), 'k': np.array([4.11149, 4.11149, 1.8329e3])},
    1250: {'Zb': np.array([1, 28.072, 1e4]), 'k': np.array([4.0128, 4.1118, 1.5973e3])}
}

# Log-log transform the data
for r, values in data.items():
    values['log_Zb'] = np.log10(values['Zb'])
    values['log_k'] = np.log10(values['k'])

# Fit piecewise linear functions (in log-log scale)
params = {}
for r, values in data.items():
    log_Zb, log_k = values['log_Zb'], values['log_k']

    # Fit two segments manually for simplicity
    # You can automate this by detecting breakpoints if necessary
    breakpoint = len(log_Zb) // 2

    slope1, intercept1 = np.polyfit(log_Zb[:breakpoint + 1], log_k[:breakpoint + 1], 1)
    slope2, intercept2 = np.polyfit(log_Zb[breakpoint:], log_k[breakpoint:], 1)

    params[r] = {
        'slope1': slope1, 'intercept1': intercept1,
        'slope2': slope2, 'intercept2': intercept2,
        'breakpoint': log_Zb[breakpoint]
    }

# Interpolation functions for slopes and intercepts
rs = np.array(list(params.keys()))
slopes1 = np.array([params[r]['slope1'] for r in rs])
intercepts1 = np.array([params[r]['intercept1'] for r in rs])
slopes2 = np.array([params[r]['slope2'] for r in rs])
intercepts2 = np.array([params[r]['intercept2'] for r in rs])
breakpoints = np.array([params[r]['breakpoint'] for r in rs])

slope1_interp = interpolate.interp1d(rs, slopes1, kind='linear', fill_value="extrapolate")
intercept1_interp = interpolate.interp1d(rs, intercepts1, kind='linear', fill_value="extrapolate")
slope2_interp = interpolate.interp1d(rs, slopes2, kind='linear', fill_value="extrapolate")
intercept2_interp = interpolate.interp1d(rs, intercepts2, kind='linear', fill_value="extrapolate")
breakpoint_interp = interpolate.interp1d(rs, breakpoints, kind='linear', fill_value="extrapolate")


# Function to calculate k for any Zb and r
def k_function(Zb, Rovert):
    if Rovert < 500:
        Rovert = 500
        print('R/t is under 500, now set to 500')
    if Rovert > 1250:
        Rovert = 1250
        print('R/t is over 1250, now set to 1250')
    Rovert = np.clip(Rovert, rs.min(), rs.max())  # Ensure r is within bounds
    log_Zb = np.log10(Zb)
    if log_Zb < breakpoint_interp(Rovert):
        log_k = slope1_interp(Rovert) * log_Zb + intercept1_interp(Rovert)
    else:
        log_k = slope2_interp(Rovert) * log_Zb + intercept2_interp(Rovert)
    return 10 ** log_k