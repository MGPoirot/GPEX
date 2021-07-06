from os.path import basename
import os
import gpxpy
import gpxpy.gpx
from scipy.stats import skewtest
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import cartopy.crs as ccrs
import cartopy.io.img_tiles as cimgt
from scipy.signal import savgol_filter

def read_gpx(path: str):
    # Read the data
    gpx_file = open(path, 'r')
    gpx = gpxpy.parse(gpx_file)
    for track in gpx.tracks:
        for segment in track.segments:
            points = [p for p in segment.points]
            speeds = []
            for i in range(1, len(points)):
                speeds.append(3.6 * points[i].speed_between(points[i - 1]))
            speeds = [speeds[0]] + speeds
    speeds = savgol_filter(speeds, 11, 1)
    for i in range(len(points)):
        points[i].speed = speeds[i]
    return points


def get_map_extent(points):
    lons = [p.longitude for p in points]
    lats = [p.latitude for p in points]
    lon_min, lon_max = min(lons), max(lons)
    lat_min, lat_max = min(lats), max(lats)
    d_lon, d_lat = lon_max - lon_min, lat_max - lat_min
    margin = 0.1

    extent = [lon_min - margin * d_lon,
             lon_max + margin * d_lon,
             lat_min - margin * d_lat,
             lat_max + margin * d_lat]

    return extent


def threshold(points):
    # Fit the data to increase normality
    result = []
    thresholds = np.linspace(6, np.median(speeds), 100)
    for t in thresholds:
        filtered_speeds = [s for s in speeds if s > t]
        _, p = skewtest(filtered_speeds)
        result.append(p)
    optim = thresholds[result.index(max(result))]
    points_in = [p for p in points if p.speed > optim]
    points_out = [p for p in points if p.speed <= optim]
    return points_in, points_out


def histplot(points_in, points_out, cmap, clim):
    cmin, cmax = clim
    speeds_in, speeds_out = [p.speed for p in points_in], [p.speed for p in points_out]
    bins = np.linspace(0, cmax * 1.5, 100)
    h = plt.hist(speeds_in, bins=bins)
    for i in range(len(bins) - 1):
        rel_col = (bins[i] - cmin) / (cmax - cmin) * 255
        rgba = cmap(int(rel_col))
        h[2].patches[i].set_color(rgba)
    plt.hist(speeds_out, bins=bins, color='r')

    plt.title(str(np.round(np.average(speeds_in), 1)) + r' km/h avg ' + '({:.0%})'.format(len(points_in) / len(points)))
    plt.xlim([0, max(bins)])
    plt.yticks([])


def scatterplot(points_in, points_out, cmap, clim):
    cmin, cmax = clim
    speeds_in, speeds_out = [p.speed for p in points_in], [p.speed for p in points_out]
    plt.scatter([p.time.timestamp() for p in points_in], speeds_in, c=speeds_in, s=1, cmap=cmap, vmin=cmin, vmax=cmax)
    plt.scatter([p.time.timestamp() for p in points_out], speeds_out, s=1, c='r')
    plt.scatter([p.time.timestamp() for p in points_in], smooth_speeds, color='k', s=1)
    plt.xlim([min([p.time.timestamp() for p in points]), max([p.time.timestamp() for p in points])])
    plt.ylim([0, cmax * 1.5])
    plt.xticks([])
    plt.ylabel('speed (km/h)')


def mapplot(points_in, points_out, cmap, clim):
    cmin, cmax = clim
    speeds_in, speeds_out = [p.speed for p in points_in], [p.speed for p in points_out]
    stamen_terrain = cimgt.Stamen('terrain-background')
    ax = plt.subplot(1, 2, 2, projection=stamen_terrain.crs)

    extent = get_map_extent(points)
    ax.set_extent(extent)

    area = (extent[1]-extent[0]) * (extent[3]-extent[2])
    detail = int(11 + np.floor(np.abs(np.log10(area / 5))))
    ax.add_image(stamen_terrain, detail)

    # plt.scatter([p.longitude for p in points_in], [p.latitude for p in points_in],
    #            c=speeds_in, s=15, transform=ccrs.Geodetic(), cmap=cmap, vmin=cmin, vmax=cmax)
    plt.scatter([p.longitude for p in points_in], [p.latitude for p in points_in],
                c=speeds_in, s=15, transform=ccrs.PlateCarree(), cmap=cmap, vmin=cmin, vmax=cmax)
    
    if len(points_out):
        # plt.scatter([p.longitude for p in points_out], [p.latitude for p in points_out],
        #            c='r', s=3, transform=ccrs.Geodetic())
        plt.scatter([p.longitude for p in points_out], [p.latitude for p in points_out],
                    c='r', s=3, transform=ccrs.PlateCarree())


fname = r'C:\Users\mgpoirot\Downloads\Lunch_Ride.gpx'
fname = r"C:\Users\mgpoirot\Downloads\Evening_Ride.gpx"
fname = r"C:\Users\mgpoirot\Downloads\Lunch_Ride(1).gpx"
fname = r"C:\Users\mgpoirot\Downloads\Morning_Ride.gpx"
fname = r"C:\Users\mgpoirot\Downloads\Night_Run.gpx"
fname = r"C:\Users\mgpoirot\Downloads\Morning_Ride(1).gpx"
#fname = r"C:\Users\mgpoirot\Downloads\Morning_Run.gpx"
#fname = r"C:\Users\mgpoirot\Downloads\Lunch_Run.gpx"
#fname = r"C:\Users\mgpoirot\Downloads\Afternoon_Ride(1).gpx"
fname = r"D:\repositories\GPEX\Afternoon_Ride.gpx"
points = read_gpx(fname)
speeds = [p.speed for p in points]

# Extract information into usable formats

points_in, points_out = threshold(points)
cmap = matplotlib.cm.get_cmap('viridis')
smooth_speeds = savgol_filter([p.speed for p in points_in], 401, 3)
clim = min(smooth_speeds), max(smooth_speeds)

# Plot
plt.figure()
plt.subplot(2, 2, 1)
histplot(points_in, points_out, cmap, clim)
plt.subplot(2, 2, 3)
scatterplot(points_in, points_out, cmap, clim)
mapplot(points_in, points_out, cmap, clim)
plt.suptitle(basename(fname))


figure_name = fname.split(os.extsep)[0] + '{}.png'
extension = ''
if os.path.isfile(figure_name.format(extension)):
    extension = 1
    while os.path.isfile(figure_name.format(str(extension))):
        extension += 1
plt.savefig(figure_name.format(str(extension)))

plt.show()
do_save = input('Do you want to save the figure? ')
if do_save not in 'Yesyes':
    os.remove(figure_name.format(str(extension)))
