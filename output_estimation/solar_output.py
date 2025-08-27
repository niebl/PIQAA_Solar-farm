import math
import os
from osgeo import gdal
import rasterio

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
    def yr_power_output(self, gti):
        # efficiency
        kWh_out = (gti * 365) * self.efficiency #get single m² efficiency
        kWh_out = kWh_out * self.size #scale to panel size

        return kWh_out

    # https://www.researchgate.net/publication/236314649_From_global_horizontal_to_global_tilted_irradiance_How_accurate_are_solar_energy_engineering_predictions_in_practice
    # based on chapter 2 of cited work
    def panel_gti(self, ghi, dni, dhi, albedo=0.15):
        # determine theta (sun incidence, by its zenith)
        solar_zenith = self.yr_avg_sun_zenith(self.lat)
        azimuth_diff = 0 #we are not assuming the solar panel to have an azimuth at this point as that would reduce efficiency
        theta = math.acos(
            math.cos(math.radians(solar_zenith)) * math.cos(math.radians(self.tilt)) +
            math.sin(math.radians(solar_zenith)) * math.sin(math.radians(self.tilt)) * math.cos(math.radians(azimuth_diff))
        ) #TODO: verify if this angle of incidence formula is correct.

        # calculate direct irradiance from the sun
        direct_tilted = dni * max(0, math.cos(theta))

        # determine transposition factors
        diffuse_transp_factor = (1 + math.cos(math.radians(self.tilt))) / 2 #this is used to transpose (tilt) DHI to the angle of the panel
        reflectance_transp_factor = (1 - math.cos(math.radians(self.tilt))) / 2 #same for the ground reflectance

        # calculate indirect irradiance
        diffuse_tilted = dhi * diffuse_transp_factor
        reflected_tilted = ghi * albedo * reflectance_transp_factor

        gti_custom = direct_tilted + diffuse_tilted + reflected_tilted
        return gti_custom

    # get the yearly average sun zenith
    # to be refined if this approximation is not sufficient
    def yr_avg_sun_zenith(self, lat):
        return lat + (23.5 / 2) #add (half of) earths tilt to the latitude
    
    def sample_irradiance(self, lat, lon, dni_name="DNI.tif", dhi_name="DIF.tif", ghi_name="GHI.tif"):
        # sample points from map
        map = arcpy.mp.ArcGISProject('CURRENT').listMaps('Map')[0]
        dni = arcpy.GetCellValue_management(map.listLayers(dhi_name)[0].dataSource, f"{lon} {lat}")
        dhi = arcpy.GetCellValue_management(map.listLayers(dni_name)[0].dataSource, f"{lon} {lat}")
        ghi = arcpy.GetCellValue_management(map.listLayers(ghi_name)[0].dataSource, f"{lon} {lat}")
        # parse results as floats
        dni = float(dni.getOutput(0))
        dhi = float(dhi.getOutput(0))
        ghi = float(ghi.getOutput(0))
        return(dni, dhi, ghi)

def main():
    panel = Solar_panel(3, 52, 7.5, efficiency=0.2, tilt=30)
    dni, dhi, ghi = panel.sample_irradiance(panel.lat, panel.lon)
    gti = panel.panel_gti(ghi, dni, dhi)
    
    panel.yr_power_output(gti)
    return

if __name__ == '__main__':
  main()