# CFos_Granule_Cells2D

**Developed for:** Beetsi

**Team:** Sélimi

**Date:** December 2022

**Software:** Fiji

### Images description

2D images taken with a x40 objective

2 channels:
  1. *CSU_405*: DAPI nuclei
  2. *CSU_642*: c-Fos cells 

With each image should be provided a .roi or .zip file containing one or multiple ROI(s).

### Plugin description

* Count nuclei based on the total volume of packed nuclei divided by as nucleus volume estimation of 250 microns^3
* Count c-Fos cells with Cellpose

### Dependencies

* **3DImageSuite** Fiji plugin

* **CLIJ** Fiji plugin

* **Cellpose* conda environment + *cyto* model

### Version history

Version 1 released on December 14, 2022.
