# Program to convert netCDF to json files
#
# Joseph B. Zambon
# 26 September 2018

# conda create --name cnaps_json python=3.5.5
# source activate cnaps_json
# conda install -c conda-forge netcdf4
# conda install scipy
# conda install matplotlib
# conda install -c conda-forge basemap
# conda install -c conda-forge proj4
# jupyter notebook --profile=nbserver

%pylab inline
from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import griddata
import matplotlib
from mpl_toolkits.basemap import Basemap

from scipy.interpolate import griddata

wrfout = '/home/ubuntu/cnaps_json/20180925.nc'
outfile = '/home/ubuntu/cnaps_json/current-wind-surface-level-gfs-0.25.json'
dx = 0.1    # Delta longitude
dy = 0.1    # Delta latitude

wrfout = Dataset(wrfout)
print(wrfout.file_format)
print(wrfout.dimensions.keys())

time = np.array(wrfout['DateTime'])
print(np.array(wrfout['DateTime']))

print(wrfout.variables.keys())

lat = np.array(wrfout['lat'])
#lat = np.squeeze(lat[0,:])
lon = np.array(wrfout['lon'])
#lon = np.squeeze(lon[:,0])
time = np.array(wrfout['DateTime'])
u_10m_tr = np.array(wrfout['u_10m_tr'])
v_10m_tr = np.array(wrfout['v_10m_tr'])

d_u_10m_tr = griddata((lon.ravel(), lat.ravel()),\
                      np.squeeze(u_10m_tr[0,:,:]).ravel(),\
                      (d_lon, d_lat), method='nearest')
d_v_10m_tr = griddata((lon.ravel(), lat.ravel()),\
                      np.squeeze(v_10m_tr[0,:,:]).ravel(),\
                      (d_lon, d_lat), method='nearest')

print(d_u_10m_tr)

figsize(25,15)
map = Basemap(projection='merc',
      resolution='l',lat_0=((np.max(lat)-np.min(lat))/2),
      lon_0=((np.max(lon)-np.min(lon))/2),
      llcrnrlon=np.min(lon),llcrnrlat=np.min(lat),
      urcrnrlon=np.max(lon),urcrnrlat=np.max(lat))
map.drawcoastlines()
map.drawcountries()
map.drawstates()
map.pcolormesh(d_lon,d_lat,d_u_10m_tr[:,:],vmin=0,vmax=10,latlon='true')

figsize(25,15)
map = Basemap(projection='merc',
      resolution='l',lat_0=((np.max(lat)-np.min(lat))/2),
      lon_0=((np.max(lon)-np.min(lon))/2),
      llcrnrlon=np.min(lon),llcrnrlat=np.min(lat),
      urcrnrlon=np.max(lon),urcrnrlat=np.max(lat))
map.drawcoastlines()
map.drawcountries()
map.drawstates()
map.pcolormesh(d_lon,d_lat,d_v_10m_tr[:,:],vmin=0,vmax=10,latlon='true')

figsize(25,15)
map = Basemap(projection='merc',
      resolution='l',lat_0=((np.max(lat)-np.min(lat))/2),
      lon_0=((np.max(lon)-np.min(lon))/2),
      llcrnrlon=np.min(lon),llcrnrlat=np.min(lat),
      urcrnrlon=np.max(lon),urcrnrlat=np.max(lat))
map.drawcoastlines()
map.drawcountries()
map.drawstates()
map.pcolormesh(d_lon,d_lat,sqrt(d_u_10m_tr[:,:]**2 +\
                                d_u_10m_tr[:,:]**2),\
               vmin=0,vmax=10,latlon='true')

d_u_10m_tr = np.flipud(d_u_10m_tr)
d_v_10m_tr = np.flipud(d_v_10m_tr)

print(time[0])

json_out = open(outfile,'w')
curr_time = str(time[0])
yyyy = curr_time[0:4]
mm = curr_time[4:6]
dd = curr_time[6:8]
hh = curr_time[8:10]
points = np.shape(d_lon[0])[0] * np.shape(d_lat[0])[0]
nx = np.shape(d_lon[0,:])[0]
ny = np.shape(d_lat[:,0])[0]

# U-component of surface wind
outstr = str("[ \n\
    {\n\
        \"header\":{\n\
            \"discipline\":0,\n\
            \"disciplineName\":\"Meteorological products\",\n\
            \"gribEdition\":2,\n\
            \"gribLength\":791507,\n\
            \"center\":8,\n\
            \"centerName\":\"NCSU CNAPS\",\n\
            \"subcenter\":0,\n\
            \"refTime\":\"" + str(yyyy) + "-" + str(mm) + "-" + str(dd) + "T" + str(hh) +":00:00.000Z\",\n\
            \"significanceOfRT\":1,\n\
            \"significanceOfRTName\":\"Start of forecast\",\n\
            \"productStatus\":0,\n\
            \"productStatusName\":\"Operational products\",\n\
            \"productType\":1,\n\
            \"productTypeName\":\"Forecast products\",\n\
            \"productDefinitionTemplate\":0,\n\
            \"productDefinitionTemplateName\":\"Analysis/forecast at horizontal level/layer at a point in time\",\n\
            \"parameterCategory\":2,\n\
            \"parameterCategoryName\":\"Momentum\",\n\
            \"parameterNumber\":2,\n\
            \"parameterNumberName\":\"U-component_of_wind\",\n\
            \"parameterUnit\":\"m.s-1\",\n\
            \"genProcessType\":2,\n\
            \"genProcessTypeName\":\"Forecast\",\n\
            \"forecastTime\":0,\n\
            \"surface1Type\":103,\n\
            \"surface1TypeName\":\"Specified height level above ground\",\n\
            \"surface1Value\":10.0,\n\
            \"surface2Type\":255,\n\
            \"surface2TypeName\":\"Missing\",\n\
            \"surface2Value\":0.0,\n\
            \"gridDefinitionTemplate\":0,\n\
            \"gridDefinitionTemplateName\":\"Latitude_Longitude\",\n\
            \"numberPoints\":" + str(points) +",\n\
            \"shape\":6,\n\
            \"shapeName\":\"Earth spherical with radius of 6,371,229.0 m\",\n\
            \"gridUnits\":\"degrees\",\n\
            \"resolution\":48,\n\
            \"winds\":\"true\",\n\
            \"scanMode\":0,\n\
            \"nx\":" + str(nx) + ",\n\
            \"ny\":" + str(ny) + ",\n\
            \"basicAngle\":0,\n\
            \"lo1\":" + str(round(min(d_lon[0,:]),2)) + ",\n\
            \"la1\":" + str(round(max(d_lat[:,0]),2)) + ",\n\
            \"lo2\":" + str(round(max(d_lon[0,:]),2)) + ",\n\
            \"la2\":" + str(round(min(d_lat[:,0]),2)) + ",\n\
            \"dx\":" + str(dx) + ",\n\
            \"dy\":" + str(dy) + ",\n\
        },\
\n")

outstr += str("        \"data\":[\n")
for j in range(0,ny):
    for i in range(0,nx):
        if j != (ny-1) or i != (nx-1):
            outstr += str("            " + str(d_u_10m_tr[j,i]) + ",\n")
        else:
            outstr += str("            " + str(d_u_10m_tr[j,i]) + "\n")

outstr += str("        ]\n\
    },\n\
")

# V-component of surface wind
            
outstr += str("    {\n\
        \"header\":{\n\
            \"discipline\":0,\n\
            \"disciplineName\":\"Meteorological products\",\n\
            \"gribEdition\":2,\n\
            \"gribLength\":791507,\n\
            \"center\":8,\n\
            \"centerName\":\"NCSU CNAPS\",\n\
            \"subcenter\":0,\n\
            \"refTime\":\"" + str(yyyy) + "-" + str(mm) + "-" + str(dd) + "T" + str(hh) +":00:00.000Z\",\n\
            \"significanceOfRT\":1,\n\
            \"significanceOfRTName\":\"Start of forecast\",\n\
            \"productStatus\":0,\n\
            \"productStatusName\":\"Operational products\",\n\
            \"productType\":1,\n\
            \"productTypeName\":\"Forecast products\",\n\
            \"productDefinitionTemplate\":0,\n\
            \"productDefinitionTemplateName\":\"Analysis/forecast at horizontal level/layer at a point in time\",\n\
            \"parameterCategory\":2,\n\
            \"parameterCategoryName\":\"Momentum\",\n\
            \"parameterNumber\":3,\n\
            \"parameterNumberName\":\"V-component_of_wind\",\n\
            \"parameterUnit\":\"m.s-1\",\n\
            \"genProcessType\":2,\n\
            \"genProcessTypeName\":\"Forecast\",\n\
            \"forecastTime\":0,\n\
            \"surface1Type\":103,\n\
            \"surface1TypeName\":\"Specified height level above ground\",\n\
            \"surface1Value\":10.0,\n\
            \"surface2Type\":255,\n\
            \"surface2TypeName\":\"Missing\",\n\
            \"surface2Value\":0.0,\n\
            \"gridDefinitionTemplate\":0,\n\
            \"gridDefinitionTemplateName\":\"Latitude_Longitude\",\n\
            \"numberPoints\":" + str(points) +",\n\
            \"shape\":6,\n\
            \"shapeName\":\"Earth spherical with radius of 6,371,229.0 m\",\n\
            \"gridUnits\":\"degrees\",\n\
            \"resolution\":48,\n\
            \"winds\":\"true\",\n\
            \"scanMode\":0,\n\
            \"nx\":" + str(nx) + ",\n\
            \"ny\":" + str(ny) + ",\n\
            \"basicAngle\":0,\n\
            \"lo1\":" + str(round(min(d_lon[0,:]),2)) + ",\n\
            \"la1\":" + str(round(max(d_lat[:,0]),2)) + ",\n\
            \"lo2\":" + str(round(max(d_lon[0,:]),2)) + ",\n\
            \"la2\":" + str(round(min(d_lat[:,0]),2)) + ",\n\
            \"dx\":" + str(dx) + ",\n\
            \"dy\":" + str(dy) + ",\n\
        },\
\n")

outstr += str("        \"data\":[\n")
for j in range(0,ny):
    for i in range(0,nx):
        if j != (ny-1) or i != (nx-1):
            outstr += str("            " + str(d_v_10m_tr[j,i]) + ",\n")
        else:
            outstr += str("            " + str(d_v_10m_tr[j,i]) + "\n")

outstr += str("        ]\n\
    }\n\
]\n\
")

json_out.write(outstr)
json_out.close()

