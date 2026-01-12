
# Engines of Continental Transformation: Africa Land Cover Dynamics (1985–2022)


This repository contains the computational framework and source code for the paper **"Engines of Continental Transformation: Accelerating Land Cover Dynamics and Systematic Transitions Across Africa (1985–2022)"**.

This project implements a **Spatio-Temporal Intensity Analysis** pipeline combining **Google Earth Engine (GEE)** for large-scale data processing and **Python** for intensity modeling and visualization.

## 📂 Repository Structure

* `gee_code/`: JavaScript algorithms for Google Earth Engine.
    * `01_matrix_generation.js`: Generates transition matrices from GLC_FCS30D data.
    * `02_spatial_mapping.js`: Maps dominant transition pathways (Deforestation, Desertification, etc.).
    * *Note: These scripts use the ROI shapefiles located in `IPCC_5_roi_shapefiles_for_africa/`.*
* `python_code/`: Python scripts for Intensity Analysis.
    * `main.py`: Performs Interval, Category, and Transition level analysis and generates figures.
    * *Usage Note: Update the `matrix_folder` variable in the script to point to your input data.*
* `Data/`: Contains sample transition matrices to test the code.
* `IPCC_5_roi_shapefiles_for_africa/`: Shapefiles used for regional stratification.
