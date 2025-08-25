import math
import os
from osgeo import gdal
import rasterio
from rasterio.plot import show
import geopandas as gpd

class Solar_panel:
    def __init__(self, size, lat, lon, efficiency=0.2, tilt=None):
        self.lat = lat #lat of panel centroid
        self.lon = lon #lon of panel centroid
        self.size = size
        self.efficiency = efficiency
        if tilt is None:
            self.tilt = 90-lat #optimum tilt if none given
        else:
            self.tilt = tilt
        return


# panel: instance of class Solar_panel
# gti: gti of solar panel adjusted to panels tilt. GTI (in kWh / m²) is daily total of yearly average.
def yr_power_output(panel, gti):
    # efficiency
    kWh_out = (gti * 365) * panel.efficiency #get single m² efficiency
    kWh_out = kWh_out * panel.size #scale to panel size

    return kWh_out

# https://www.researchgate.net/publication/236314649_From_global_horizontal_to_global_tilted_irradiance_How_accurate_are_solar_energy_engineering_predictions_in_practice
# based on chapter 2 of cited work
def custom_gti(panel, ghi, dni, dhi, albedo=0.15):
    # determine theta (sun incidence, by its zenith)
    solar_zenith = yr_avg_sun_zenith(panel.lat)
    azimuth_diff = 0 #we are not assuming the solar panel to have an azimuth at this point as that would reduce efficiency
    theta = math.acos(
        math.cos(math.radians(solar_zenith)) * math.cos(math.radians(panel.tilt)) +
        math.sin(math.radians(solar_zenith)) * math.sin(math.radians(panel.tilt)) * math.cos(math.radians(azimuth_diff))
    ) #TODO: verify if this angle of incidence formula is correct.

    # calculate direct irradiance from the sun
    direct_tilted = dni * max(0, math.cos(theta))

    # determine transposition factors
    diffuse_transp_factor = (1 + math.cos(math.radians(panel.tilt))) / 2 #this is used to transpose (tilt) DHI to the angle of the panel
    reflectance_transp_factor = (1 - math.cos(math.radians(panel.tilt))) / 2 #same for the ground reflectance

    # calculate indirect irradiance
    diffuse_tilted = dhi * diffuse_transp_factor
    reflected_tilted = ghi * albedo * reflectance_transp_factor

    gti_custom = direct_tilted + diffuse_tilted + reflected_tilted
    return gti_custom

def sample_gtif(file, lat, lon):
    raster = rasterio.open(file)
    sample = raster.sample([(lon,lat)])
    value = next(sample)
    return value[0]

def sample_irradiance(datapath, lat, lon):
    dni = sample_gtif(os.path.join(datapath, "DNI.tif"), lat, lon)
    dhi = sample_gtif(os.path.join(datapath, "DIF.tif"), lat, lon)
    ghi = sample_gtif(os.path.join(datapath, "GHI.tif"), lat, lon)
    return(dni, dhi, ghi)

# get the yearly average sun zenith
# to be refined if this approximation is not sufficient
def yr_avg_sun_zenith(lat):
    return lat + (23.5 / 2) #add (half of) earths tilt to the latitude

def main():
    panel = Solar_panel(3, 52, 7.5, efficiency=0.2, tilt=30)

    dni, dhi, ghi = sample_irradiance(".", panel.lat, panel.lon)
    print(dni)
    print(ghi)
    print(dhi)

    print(sample_gtif("GTI.tif", panel.lat, panel.lon))

    gti = custom_gti(panel, ghi, dni, dhi)
    print(gti)
    print(yr_power_output(panel,gti))
    return

if __name__ == '__main__':
  main()