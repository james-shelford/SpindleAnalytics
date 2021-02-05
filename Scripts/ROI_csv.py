"""
save_rois_csv.py 
created by: Erick Martins Ratamero
date: 26/09/18
From an open image in Fiji with ROIs,
save each individual ROI as a CSV file.
"""

from ij.plugin.frame import RoiManager
import os


# returns an instance of ROI Manager (creates one if 
# it doesn't exist)
def get_roi_manager(new=False):
    rm = RoiManager.getInstance()
    if not rm:
        rm = RoiManager()
    if new:
        rm.runCommand("Reset")
    return rm



# getting the ROIs generated
rm = get_roi_manager()

rois = rm.getRoisAsArray()
for roi in rois:
	polygon = roi.getPolygon()
	xp = polygon.xpoints
	yp = polygon.ypoints
	name = roi.getName()+'.csv'
	f = open(name, 'w')
	for counter in range(len(xp)):
		f.write(','.join([str(xp[counter]), str(yp[counter])])+'\n')
	f.close()