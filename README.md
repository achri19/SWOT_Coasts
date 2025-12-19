# *SWOT_Coasts*

### Several workflows to extract SWOT data for evaluation against in-situ data in coastal areas


# Requirements:
## Python packages
```
conda env create -f requirements.yml

conda activate swot
```

## Configuration file
A CSV of AOIs and bounding boxes will be read into the first two workflows below
```
examples_template.csv 
```
# Workflows

## 1. Download SWOT data (uses earthaccess)
   Answer prompts to determine which products to download and provide a bounding box.
   ```
   python SWOT_download_files.py --interactive
   ```

   The examples_template.csv file will be read and you will be asked to choose which AOI you want to use. The corresponding bounding box will be used to search for SWOT data. 
   
   You can choose from LR, HR Raster, HR PIXC, HR RiverSP, HR Lake SP, and HR PIXCVec and Version C or D. However, only the first 3 are included in the other workflows of this repo (as of Dec 19, 2025). 
   
   If you choose LR, you will be given the opportunity to download L3 LR products through AVISO (https://www.aviso.altimetry.fr/en/data/products/sea-surface-height-products/global/swot-l3-ocean-products.html), but you will need an account and this step will be time consuming.

   
## 2. Find track/frame information for a set of Lat/Lon coordinates
   - Determine the SWOT pass, scene, and tiles that overlap a given point or bounding box
   - Provide path to either: 
      - CSV with point coordinates with columns: name, latitude, and longitude 
      - CSV with bounding box coordinates with columns: aoi, minx, miny, maxx, maxy
   
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
