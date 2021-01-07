import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy import special
import matplotlib.pyplot as plt
from scipy import integrate

def cumgauss(x, mu, sigma):
    return 0.5 * (1 + special.erf((x-mu)/(np.sqrt(2)*sigma)))
	
e = 11567.3
estartvar = 11500  # 11120 
eendvar = 11577.1 # 11231.7

width = 5.31 #5.25 for L3
	
folder = r'Z:\Standards\PtL\PtO2\output_Pt_L3_PtO2\Export'
file_raw = r'PtO2_sam_matrix_Both_normalised_all.dat'
	
import_data_raw = pd.read_csv(folder+'/'+file_raw, sep='\t', header=0)

energy = import_data_raw['E']
mu = import_data_raw['0']

cumgauss_step = cumgauss(energy, e, width)

plt.plot(energy, mu)
plt.plot(energy, cumgauss_step)

diff = mu - cumgauss_step

diff_df = pd.DataFrame()
diff_df['E'] = energy
diff_df['diff'] = diff

diff_df = diff_df[(diff_df['E'] >= float(estartvar)) & (diff_df['E'] <= float(eendvar))]
y_int = integrate.cumtrapz(diff_df['diff'], diff_df['E'], initial=0)

plt.plot(diff_df['E'], diff_df['diff'])
plt.show()

plt.plot(y_int)
plt.show()

print('integral :', y_int[-1])