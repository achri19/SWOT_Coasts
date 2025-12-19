# *SWOT_Coasts*

### Several workflows to extract SWOT data for evaluation against in-situ data in coastal areas


# Requirements:
## Python packages
```
conda env create -f requirements.yml

conda activate swot
```

## Configuration files
CSVs of AOIs and bounding boxes will be read into the first two workflows below:

   - *aoi_template.csv* contains rows for areas of interest: aoi, minx, miny, maxx, maxy
   - *point_template.csv* contains rows with individual points: name, aoi, latitude, longitude

A demo location (Amerada Pass in the Atchafalaya Riva, Louisiana) is included. Users are encouraged to add their AOIs and individual points to each CSV to download and process SWOT data for their study areas


# Workflows

## 1. Download SWOT data (uses earthaccess)
   ```
   python SWOT_download_files.py --interactive
   ```

   - The *aoi_template.csv* file will be read and you will be asked to choose which AOI you want to use. The corresponding bounding box will be used to search for SWOT data. 
   
   - You can choose from LR, HR Raster, HR PIXC, HR RiverSP, HR Lake SP, and HR PIXCVec and Version C or D. However, only the first 3 products are included in the other workflows of this repo (as of Dec 19, 2025). 
   
   - If you choose LR, you will be given the opportunity to download L3 LR products through AVISO (https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/global/swot-l3-ocean-products.html), but you will need an account and this step will be time consuming.

   
## 2. Find track/frame information for a set of Lat/Lon coordinates
   - Determine the SWOT pass, scene, and tiles that overlap your AOI: both the entire bounding box (given in aoi_template.csv) and indivudal points (point_template.csv). Determining the passes, scenes, and tiles reduces the processing time of later workflows.
   - This workflow is time consuming, but only needs to be run for each AOI and list of points one time.
   
   ```
   python SWOT_get-pass-scene-tile.py
   ```

## Extract SWOT HR Raster data around a point
   - Extract SWOT data within a given search radius (km)
   - Updated water surface elevation (WSE) values using a tide-free EGM08 geoid (SWOT geoid is mean tide EGM08 geoid)
   - Add "good" flag that indicates which pixels have wse_qual == 0 and only with 10-60km on either side of the nadir.
   - If a watermask is available, also flags pixels are being within or not within the water mask
   
   ```
   python SWOT_L2_HR_Raster_ExtractNearPoint --interactive
   ```
