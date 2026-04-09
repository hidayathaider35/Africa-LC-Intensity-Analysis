# Engines of Continental Transformation: Africa Land Cover Dynamics (1985–2022)

This repository contains the computational framework and source code for the paper **"Significant Land Cover Transitions and Regional Acceleration at the Continental Scale of Africa Over the Last Four Decades"**.

This project implements a **Spatio-Temporal Intensity Analysis** pipeline combining **Google Earth Engine (GEE)** for large-scale data processing and **Python** for intensity modeling and visualization.

## 📂 Repository Structure

* `gee_code/`: JavaScript algorithms for Google Earth Engine.
    * `01_matrix_generation.js`: Generates transition matrices from GLC_FCS30D data.
    * `02_spatial_mapping.js`: Maps dominant transition pathways (Deforestation, Desertification, etc.).
    * *Note: These scripts use the ROI shapefiles located in `IPCC_5_roi_shapefiles_for_africa/`.*

* `python_code/`: Python scripts for Intensity Analysis.
    * `main.py`: Performs Interval, Category, and Transition level analysis with regional decomposition. Generates all main-text results (Figures 2–6) in .csv and .png.
    * `supplementary_generator.py`: Generates per-region supplementary figures (Figures S1–S5), including interval-, category-, and transition-level analyses for each IPCC sub-region, and temporal trajectory plots.
    * *Usage Note: Update the `matrix_folder` and `output_folder` variables in each script to point to your local directories.*

* `Data/`: Contains sample transition matrices to test the code.

* `IPCC_5_roi_shapefiles_for_africa/`: Shapefiles for the five IPCC AR5 African climate reference regions used for regional stratification.

## 🔧 Requirements

* **Google Earth Engine**: A GEE account is required to run the JavaScript scripts.
* **Python 3.8+** with the following packages:
    ```
    numpy
    pandas
    matplotlib
    seaborn
    ```

## 🚀 Usage

### Step 1: Generate Transition Matrices (GEE)
Run `gee_code/01_matrix_generation.js` in the GEE Code Editor. This produces CSV transition matrices for each region and time interval, exported to Google Drive.

### Step 2: Main Analysis (Python)
```bash
python python_code/main.py
```
Update `matrix_folder` to point to the downloaded CSV files. Outputs: main-text figures and CSV tables for interval, category, and transition levels.

### Step 3: Supplementary Figures (Python)
```bash
python python_code/supplementary_generator.py
```
Update `matrix_folder` and `output_folder` paths. Outputs: per-region intensity analysis figures and temporal trajectory plots.

## 📊 Dataset

This study uses the **GLC_FCS30D** dataset (30 m, 1985–2022):
* Zhang et al. (2024). *Earth System Science Data*, 16, 1353–1381. [doi:10.5194/essd-16-1353-2024](https://doi.org/10.5194/essd-16-1353-2024)
* GEE catalog: `projects/sat-io/open-datasets/GLC-FCS30D`

## 📄 Citation

If you use this code, please cite:

```[Ullah, H., Kalisa, W., Ali, S., Kong, D., Zhang, J., 2026. Significant Land Cover Transitions and Regional Acceleration at the Continental Scale of Africa over the Last Four Decades. Sensors. 26, 2318. https://doi.org/10.3390/S26082318```
