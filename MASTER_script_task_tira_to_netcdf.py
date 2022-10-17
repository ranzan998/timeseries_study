import os

import pandas as pd
import xarray as xr
import numpy as np

os.system('rm -rf  tide.out tide1.out tide2.out tide.in')
os.system('rm -rf tira.ctl tira.o format.o ')

# to list all the files in this directory
file_name = sorted(os.listdir('/home/scilab/Desktop/shankar_fifth/1datasets/first_dsets/'))

print(file_name)
#to count the length of each of the above files, to be used by 'nl' and NREC of tira1m.ctl file later
len_file = []
for f in file_name:
    data = xr.open_dataset('/home/scilab/Desktop/shankar_fifth/1datasets/first_dsets/'+f)
    uvel = data["SPEED"]
    nl = len(uvel)
    len_file.append(nl)
len_file.append(8895)
print(len_file)
####################################################
#change the NREC value in tira1m.ctl for the next run
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
#to rewrite tira1m.ctl to tira.ctl as a new file
def copy(old_name, new_name):
    with open(old_name) as file:
        s = file.read()
    with open(new_name, 'w') as file:
        file.write(s)
        file.close()
###############################
#  X    X       X       X       X       X main body X    X       X       X       X       X#
c = 0
for f in file_name:
    data = xr.open_dataset('/home/scilab/Desktop/shankar_fifth/1datasets/first_dsets/' + f)
    uvel1 = data["SPEED"] * np.sin(np.deg2rad(data["DIRECTION"]))
    vvel1 = data["SPEED"] * np.cos(np.deg2rad(data["DIRECTION"]))

    uvel = [(uvel1[i] if -100<uvel1[i]<100 else 0) for i in range(len(uvel1))]
    vvel = [(vvel1[i] if -100<vvel1[i]<100 else 0) for i in range(len(uvel1))]

    nl = len(uvel)
    copy('tira1m.ctl', 'tira.ctl')

    date_numpy = data["TAXIS"].dt.strftime("%d-%^b-%Y %H:%M:%S").values

    #write the tide.txt files
    with open('tide.txt', 'w') as file:
        for i in range(len(uvel)):
            if len(uvel) < 10000:
                file.write(' ' + str(date_numpy[i]) + ' /' + str(i + 1).rjust(5, " ") + ': ' + str(
                    np.round(float(uvel[i]), 2)).rjust(6) + "\n")
            else:
                file.write(' ' + str(date_numpy[i]) + ' /' + str(i + 1).rjust(6, " ") + ': ' + str(
                    np.round(float(uvel[i]), 2)).rjust(6) + "\n")
    file.close()

    #run TASK for zonal current
    os.system("gfortran format.f -o format.o")
    os.system("gfortran tira.f -o tira.o")
    os.system("./format.o; ./tira.o")

    #write the tira.pri file > convert to csv using sed
    with open("tira.pri", 'r') as file:
        line_count = len(file.readlines())  #count the line to delete the exact unnecessary lines

    os.system(
        "sed -n '%s,%s p' tira.pri > tira1.txt; sed '2,3d' tira1.txt > tira11.txt;sed 's/  */,/g' tira11.txt > tira_uvel.csv" % (
        str(line_count - 26), str(line_count)))
    df = pd.read_csv('tira_uvel.csv')
    df1= df.sort_values('H',ascending=False).head(10)
    df1.to_csv(f+'tira_uvel.csv')

    #remove the files
    os.system("cat tide.out >> tide1.out; rm -rf tira1.txt tira11.txt tira_uvel.csv")


    #repeating steps for meridional current
    with open('tide.txt', 'w') as file:
        for i in range(len(uvel)):
            if len(uvel) < 10000:
                file.write(' ' + str(date_numpy[i]) + ' /' + str(i + 1).rjust(5, " ") + ': ' + str(
                    np.round(float(vvel[i]), 2)).rjust(6) + "\n")
            else:
                file.write(' ' + str(date_numpy[i]) + ' /' + str(i + 1).rjust(6, " ") + ': ' + str(
                    np.round(float(vvel[i]), 2)).rjust(6) + "\n")
    file.close()
    os.system("gfortran format.f -o format.o")
    os.system("gfortran tira.f -o tira.o")
    os.system("./format.o; ./tira.o")

    os.system(
        "sed -n '%s,%s p' tira.pri > tira2.txt; sed '2,3d' tira2.txt > tira22.txt;sed 's/  */,/g' tira22.txt > tira_vvel.csv" % (
        str(line_count - 26), str(line_count)))
    df = pd.read_csv('tira_vvel.csv')
    df1 = df.sort_values('H',ascending=False).head(10)
    df1.to_csv(f+'tira_vvel.csv')
    os.system("cat tide.out >> tide2.out; rm -rf tira2.txt tira22.txt tira_vvel.csv")


    U_RAW =  np.array(np.round(np.loadtxt("tide1.out")[:, 5].astype(np.float32), 2))

    U_RES = np.loadtxt("tide1.out")[:, 9]
    U_TIDE = np.loadtxt("tide1.out")[:, 8]
    V_RAW = np.loadtxt("tide2.out")[:, 5]
    V_RES = np.loadtxt("tide2.out")[:, 9]
    V_TIDE = np.loadtxt("tide2.out")[:, 8]

    depth =[1]
    dep = xr.Dataset({
        'depth': xr.DataArray(
            data=depth,
            dims=['depth'],
            attrs={'long_name': 'Depth (m)', 'units': 'meters', 'positive': 'down', 'point_spacing': 'even',
                   'axis': 'Z'}
        )
    })
    lat = xr.DataArray([9.1931])
    lon =xr.DataArray([79.3104])
    ds = xr.Dataset({
        'URAW': xr.DataArray(
            data=np.reshape(U_RAW,(len(U_RAW),1,1,1)),  # enter data here
            dims=['time','depth','lon','lat'],
            coords={'time': data["TAXIS"].values,'depth':dep["depth"].values,'lon':lon.values,'lat':lat.values},
            attrs={'long_name': 'u raw', 'units': 'cm/s'}
        ),
        'VRAW': xr.DataArray(
            data=np.reshape(V_RAW,(len(V_RAW),1,1,1)),  # enter data here
            dims=['time','depth','lon','lat'],
            coords={'time': data["TAXIS"].values,'depth':dep["depth"].values,'lon':lon.values,'lat':lat.values},
            attrs={'long_name': 'v raw', 'units': 'cm/s'}
        ),
        'UTIDE': xr.DataArray(
            data=np.reshape(U_TIDE,(len(U_TIDE),1,1,1)),  # enter data here
            dims=['time','depth','lon','lat'],
            coords={'time': data["TAXIS"].values,'depth':dep["depth"].values,'lon':lon.values,'lat':lat.values},
            attrs={'long_name': 'u tide', 'units': 'cm/s'}
        ),
        'VTIDE': xr.DataArray(
            data=np.reshape(V_TIDE,(len(V_TIDE),1,1,1)),  # enter data here
            dims=['time','depth','lon','lat'],
            coords={'time': data["TAXIS"].values,'depth':dep["depth"].values,'lon':lon.values,'lat':lat.values},
            attrs={'long_name': 'v tide', 'units': 'cm/s'}
        ),
        'URES': xr.DataArray(
            data=np.reshape(U_RES,(len(U_RES),1,1,1)),  # enter data here
            dims=['time','depth','lon','lat'],
            coords={'time': data["TAXIS"].values,'depth':dep["depth"].values,'lon':lon.values,'lat':lat.values},
            attrs={'long_name': 'u residue', 'units': 'cm/s'}
        ),
        'VRES': xr.DataArray(
            data=np.reshape(V_RES,(len(U_RES),1,1,1)),  # enter data here
            dims=['time','depth','lon','lat'],
            coords={'time': data["TAXIS"].values,'depth':dep["depth"].values,'lon':lon.values,'lat':lat.values},
            attrs={'long_name': 'v residue', 'units': 'cm/s'})

    },
        attrs={'author': 'ranjan'}
    )

    ds.depth.attrs['long_name'] = 'Depth (m)'
    ds.depth.attrs['units'] = 'meters'
    ds.depth.attrs['positive'] = 'down'
    ds.depth.attrs['point_spacing'] = 'even'
    ds.depth.attrs['axis'] = 'Z'
    
    ds.time.attrs['long_name'] = 'time'
    ds.time.attrs['time_origin'] = "1901-01-15"
    ds.time.attrs['axis'] = 'T'
    ds.time.encoding['units'] = 'minutes since 1901-01-15'

    ds = ds.convert_calendar("standard", dim='time', use_cftime=True)
    
    ds.to_netcdf('/home/scilab/Desktop/shankar_fifth/2detided_all/detided_'+f, format='NETCDF4')
    os.system("rm -rf tide1.out tide2.out")

    if len_file[c + 1] < 10000 and len_file[c] < 10000:
        inplace_change('/home/scilab/Desktop/shankar_fifth/2detided_all/tira1m.ctl', str(nl ), str(len_file[c + 1]))

    elif len_file[c + 1] < 10000 and len_file[c] > 10000:
        inplace_change('/home/scilab/Desktop/shankar_fifth/2detided_all/tira1m.ctl', str(nl ), ' '+str(len_file[c + 1]))

    elif len_file[c + 1] > 10000 and len_file[c] < 10000:
        inplace_change('/home/scilab/Desktop/shankar_fifth/2detided_all/tira1m.ctl', ' '+str(nl ), str(len_file[c + 1]))

    elif len_file[c + 1] > 10000 and len_file[c] > 10000:
        inplace_change('/home/scilab/Desktop/shankar_fifth/2detided_all/tira1m.ctl', str(nl ), str(len_file[c + 1]))
    else:
        pass

    c += 1
    print(c)
