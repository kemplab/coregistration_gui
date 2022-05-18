## Coregistration

For assistance, contact Tanya Roysam _roysamt@gatech.edu_

### Overview

Graphical User Interface ('GUI') for aligning and analyzing confocal images and matrix-assisted laser 
desorption/ionization mass spectrometry images ('MALDI' MSI). 

### Features

* Parse ImzML files
* Align confocal images with MALDI mass spectrometry ion data
* Segment confocal images and overlay with aligned MALDI ion data 
* Calculate extensive cell network metrics

### Usage
All six .py files must be in the same directory.

Run ```coregistration.py``` by double-clicking or use the following bash command:
```bash
python coregistration.py
```
**Coregistration requires a .csv file of raw MALDI data in a specific format; run ImzML parsing beforehand and do not 
edit this file.**

* Navigate between pages with the 'Previous Page' and 'Next Page' buttons
* Progress between steps and following inputs with the 'Continue' button
* Input values and navigate to files using interactive file dialogs
* Explore images by zooming in, panning around, and saving snapshots using the toolbar beneath image outputs.
* Keep a record of your work and troubleshoot using session logs saved in the output directory under a timestamp
* Store all outputs in organized project directory

**Tip**: If you click the 'Continue' button and nothing happens, check that you clicked 'Save' in all applicable fields
and entered valid inputs in all required fields.


### Dependencies

* Python >= 3.8
* opencv-python >= 3.0
* numpy
* matplotlib
* scipy
* scikit-image
* sklearn
* pyimzml
* pandas
* pillow _(**not** PIL)_