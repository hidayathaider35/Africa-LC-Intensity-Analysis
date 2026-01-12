
# Engines of Continental Transformation: Africa Land Cover Dynamics (1985–2022)


This repository contains the computational framework and source code for the paper **"Engines of Continental Transformation: Accelerating Land Cover Dynamics and Systematic Transitions Across Africa (1985–2022)"**.

This project implements a **Spatio-Temporal Intensity Analysis** pipeline combining **Google Earth Engine (GEE)** for large-scale data processing and **Python** for intensity modeling and visualization.

## 📂 Repository Structure
* 
*   `gee_code/`: JavaScript algorithms for Google Earth Engine.
    *   `01_matrix_generation.js`: Generates transition matrices from GLC_FCS30D data.
    *   `02_spatial_mapping.js`: Maps dominant transition pathways (Deforestation, Desertification, etc.).
         (notes: please use the "IPCC_5_roi_shapefiles_for_africa files/" as a roi)
*   `python_code/`: Python scripts for Intensity Analysis.
    *   `main.py`: Performs Interval, Category, and Transition level analysis and generates figures.
	         (notes: for the variable “matrix_folder”, please change to  the folder where all the matrices are downloaded using the 01_matrix_generation.js code）		
*   `data/`: Contains sample transition matrices for the 5 IPCC African sub-regions.

## 🚀 Getting Started

### Prerequisites
*   A Google Earth Engine account.
*   Python 3.8+
*   Required Python libraries: `pandas`, `numpy`, `matplotlib`, `seaborn`

### Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/hidayathaider35/Africa-LC-Intensity-Analysis.git
   cd Africa-LC-Intensity-Analysis
