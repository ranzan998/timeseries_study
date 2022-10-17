import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os

file_path = sorted(os.listdir('/home/scilab/Desktop/shankar_fifth/4residuals/'))
file_names = [f for f in file_path if f.endswith('.nc')]
file_names.append(file_names[0])
def inplace_change(filename, old_string, new_string):
    with open(filename) as file:
        s = file.read()
        if old_string not in s:
            print('"{old_string}" not found in {filename}.'.format(**locals()))
            return

    with open(filename, 'w') as file:
        print('Changing "{old_string}" to "{new_string}" in {filename}'.format(**locals()))
        s = s.replace(old_string, new_string)
        file.write(s)

len_file = []
for f in file_names:
    data = xr.open_dataset('/home/scilab/Desktop/shankar_fifth/4residuals/'+f)
    uvel = data["URES1"]
    nl = len(uvel)
    len_file.append(nl)
len_file.append(8895)
print(len_file)
print(file_names)
c = 0
for f in file_names:
    print(f[-7:-3])

for f in file_names:
    os.system('/home/scilab/Desktop/shankar_fifth/4residuals/rotary_run.sh')
    os.system('mv ' + f[-7:-3]+'_ureal.nc ' + f[-7:-3]+'_uimag.nc ' + f[-7:-3]+'_vreal.nc ' + f[-7:-3]+'_vimag.nc /home/scilab/Desktop/shankar_fifth/4residuals/test_wave')
    os.system('conda run -n FERRET ferret -script wavelet_cw_ccw.jnl')
    os.system('mv '+ f[-7:-3]+'_nc_residuals_wavelet_cwccw.nc /home/scilab/Desktop/shankar_fifth/5cwccw')
    os.system('rm ./test_wave/' + f[-7:-3]+'_ureal.nc ./test_wave/' + f[-7:-3]+'_uimag.nc ./test_wave/' + f[-7:-3]+'_vreal.nc ./test_wave/' + f[-7:-3]+'_vimag.nc')
    os.system('rm wavetest_uimag.o wavetest_ureal.o wavetest_vimag.o wavetest_vreal.o')

    inplace_change('/home/scilab/Desktop/shankar_fifth/4residuals/wavetest_ureal.f', str(len_file[c]),str(len_file[c+1]))
    inplace_change('/home/scilab/Desktop/shankar_fifth/4residuals/wavetest_uimag.f', str(len_file[c]), str(len_file[c+1]))
    inplace_change('/home/scilab/Desktop/shankar_fifth/4residuals/wavetest_vreal.f', str(len_file[c]), str(len_file[c+1]))
    inplace_change('/home/scilab/Desktop/shankar_fifth/4residuals/wavetest_vimag.f', str(len_file[c]), str(len_file[c+1]))

    inplace_change('/home/scilab/Desktop/shankar_fifth/4residuals/wavetest_ureal.f',f[-7:-3],file_names[c+1][-7:-3])
    inplace_change('/home/scilab/Desktop/shankar_fifth/4residuals/wavetest_uimag.f', f[-7:-3],
                   file_names[c + 1][-7:-3])
    inplace_change('/home/scilab/Desktop/shankar_fifth/4residuals/wavetest_vreal.f', f[-7:-3],
                   file_names[c + 1][-7:-3])
    inplace_change('/home/scilab/Desktop/shankar_fifth/4residuals/wavetest_vimag.f', f[-7:-3],
                   file_names[c + 1][-7:-3])

    inplace_change('/home/scilab/Desktop/shankar_fifth/4residuals/wavelet_cw_ccw.jnl',f[-7:-3],file_names[c+1][-7:-3])

    c += 1



