import math
import os
from osgeo import gdal
import numpy as np
import arcpy

class Solar_panel:
    def __init__(self, size, lat, lon, cost, efficiency=0.2, tilt=None):
        self.lat = lat #lat of panel centroid
        self.lon = lon #lon of panel centroid
        self.size = size
        self.efficiency = efficiency
        self.dni = None
        self.dhi = None
        self.ghi = None
        self.cost = cost
        if tilt is None:
            self.tilt = 90-lat #optimum tilt if none given
        else:
            self.tilt = tilt
        return
    
    # panel: instance of class Solar_panel
    # gti: gti of solar panel adjusted to panels tilt. GTI (in kWh / m²) is daily total of yearly average.
    def yr_power_output(self, gti):
        # efficiency
        out = (gti * 365) * self.efficiency #get single m² efficiency
        out = out * self.size #scale to panel size

        return out

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
    
    def set_irradiance_maps(self, dni, dhi, ghi):
        # set irradiance maps
        self.dni_map = dni
        self.dhi_map = dhi
        self.ghi_map = ghi

    def sample_irradiance(self, lat, lon, cache=None, dni_name="DNI.tif", dhi_name="DIF.tif", ghi_name="GHI.tif"):
        # this method is too slow for large number of panels due to geoprocessing overhead
        # keeping it around for reference.
        if self.dni_map and self.dhi_map and self.ghi_map:
            # if irradiance maps already exist, do the more efficient step
            # this is still rather slow
            location = f"{lon} {lat}"
            dni = arcpy.GetCellValue_management(self.dni_map, location)
            dhi = arcpy.GetCellValue_management(self.dhi_map, location)
            ghi = arcpy.GetCellValue_management(self.ghi_map, location)
            # parse results as floats
            dni = float(dni.getOutput(0))
            dhi = float(dhi.getOutput(0))
            ghi = float(ghi.getOutput(0))
            return(dni, dhi, ghi)
        else:
            # else, do the old method. this is very slow. find a different way than re-listing the layers each time
            # sample points from map
            arcgis_map = arcpy.mp.ArcGISProject('CURRENT').listMaps('Map')[0]
            dni = arcpy.GetCellValue_management(arcgis_map.listLayers(dni_name)[0].dataSource, f"{lon} {lat}")
            dhi = arcpy.GetCellValue_management(arcgis_map.listLayers(dhi_name)[0].dataSource, f"{lon} {lat}")
            ghi = arcpy.GetCellValue_management(arcgis_map.listLayers(ghi_name)[0].dataSource, f"{lon} {lat}")
            # parse results as floats
            dni = float(dni.getOutput(0))
            dhi = float(dhi.getOutput(0))
            ghi = float(ghi.getOutput(0))
            return(dni, dhi, ghi)


def script_tool(input_panels, DNI_src, DHI_src, GHI_src, kwh_price, DNI_unit="kWh", DHI_unit="kWh", GHI_unit="kWh"):
    arcpy.SetProgressor("default", "beginning solar power output calculation")
    
    ## Loop through all panels in input feature group and create panel objects
    farm_panels = []
    keys = ["Shape@","Shape_Length","PANEL_COST","PANEL_EFF","TILT_DEG","AZIMUTH_DEG","PANEL_AREA_M2"]

    arcpy.SetProgressor("default", "progressing solar panel features")
    with arcpy.da.SearchCursor(input_panels, keys) as cursor:
        for row in cursor:
            geom = row[0].projectAs(arcpy.SpatialReference(4326)) #reproject in case the previous toolbox used 3857 or something
            centroid = geom.centroid
            new_panel = Solar_panel(
                size=row[6],
                lat=centroid.Y,lon=centroid.X,
                efficiency=row[3],
                tilt=row[4],
                cost=row[2]
                )
            farm_panels.append(new_panel)

    ## sample irradiances for power output calculation
    # (this can't be done by class methods because it's slow one by one)

    arcpy.SetProgressor("default", "sampling irradiance values for each solar panel")
    # Prepare rasters
    def get_raster_array(irr_map):
        raster = arcpy.Raster(irr_map)
        return raster, arcpy.RasterToNumPyArray(raster)
    
    dni_raster, dni_array = get_raster_array(DNI_src)
    dhi_raster, dhi_array = get_raster_array(DHI_src)
    ghi_raster, ghi_array = get_raster_array(GHI_src)

    # For each point, convert x,y to row, col
    def transform_point_coords(x, y, raster):
        xMin = raster.extent.XMin
        yMax = raster.extent.YMax
        cellWidth = raster.meanCellWidth
        cellHeight = raster.meanCellHeight
        col = int((x - xMin) / cellWidth)
        row = int((yMax - y) / cellHeight)
        return row, col

    for panel in farm_panels:
        # sample irradiance values and set them in the object again
        row, col = transform_point_coords(panel.lon, panel.lat, dni_raster)
        panel.dni = dni_array[row,col]
        row, col = transform_point_coords(panel.lon, panel.lat, dhi_raster)
        panel.dhi = dhi_array[row,col]
        row, col = transform_point_coords(panel.lon, panel.lat, ghi_raster)
        panel.ghi = ghi_array[row,col]
        
        
    ## Access functions to calculate power output of those
    arcpy.SetProgressor("default", "calculating power output of all solar panels")
    total_yr_kwh = 0
    for panel in farm_panels:
        # this samples irradiance values and calculates the global irradiance for the panel at its given tilt
        gti = panel.panel_gti(panel.ghi, panel.dni, panel.dhi)
        panel_yr_kwh = panel.yr_power_output(gti)
        total_yr_kwh = panel_yr_kwh + total_yr_kwh

    ## Access functions to calculate costs of solar panels
    total_farm_cost = 0
    for panel in farm_panels:
        total_farm_cost = total_farm_cost + panel.cost
    total_farm_cost_string = f"total cost of solar panels: {total_farm_cost}€"

    ## Calculate time-to-payoff for solar panel costs
    # here a year is exactly 365 days, as that is used in the total_yr_kwh method
    yrs_for_payoff = total_farm_cost / ((total_yr_kwh)*kwh_price)
    days_for_payoff = yrs_for_payoff * 365
    yrs_flr = math.floor(yrs_for_payoff)
    days_remainder = math.ceil(days_for_payoff - (yrs_flr * 365))
    payoff_string = f"time to payoff: {yrs_flr} years, {days_remainder} days"

    ## Return values to user
    yr_kwh = round(total_yr_kwh, 2)
    kwh_out_string = f"total yearly power output: {yr_kwh} kWh"
    arcpy.AddMessage(kwh_out_string)
    arcpy.AddMessage(total_farm_cost_string)
    arcpy.AddMessage(payoff_string)

    return

if __name__ == "__main__":
    # get parameters

    # so far, most of the relevant infos are carried in the attribute table
    param_input_panels = arcpy.GetParameterAsText(0) 
   
    param_DNI = arcpy.GetParameterAsText(1) # tif input of irradiance maps
    #param_DNI_unit = arcpy.GetParameter(2) # whether it's Wh or kWh. default to kWh
    param_DHI = arcpy.GetParameterAsText(2)
    #param_DHI_unit = arcpy.GetParameter(4)
    param_GHI = arcpy.GetParameterAsText(3)
    #param_GHI_unit = arcpy.GetParameter(6)
    param_kwh_price = arcpy.GetParameterAsText(4)

    param_kwh_price = float(param_kwh_price)

    #TODO add support for files on disk

    script_tool(
        input_panels=param_input_panels, 
        DNI_src=param_DNI, #param_DNI_unit, 
        DHI_src=param_DHI, #param_DHI_unit, 
        GHI_src=param_GHI, #param_GHI_unit,
        kwh_price=param_kwh_price)
    #arcpy.SetParameterAsText(2, "Result")
