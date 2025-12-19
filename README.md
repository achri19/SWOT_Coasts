# *SWOT_Coasts*

### Several workflows to extract SWOT data for evaluation against in-situ data in coastal areas


# Requirements:
## Python packages
```
conda env create -f requirements.yml

conda activate swot
```

## Configuration file
```
examples_template.csv provides an example and will be read into the first two workflows below
```
# Workflows

## 1. Download SWOT data (uses earthaccess)
   Answer prompts to determine which products to download and provide a bounding box.
   ```
   python SWOT_download_files.py --interactive
   ```

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
