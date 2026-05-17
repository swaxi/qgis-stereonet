# qgis-stereonet v0.3.02
 WAXI QFIELD Fork of steronet plugin

# Source
 All the hard work was carried out by Joe Kington: https://github.com/joferkington/mplstereonet and Daniel Childs: https://github.com/childsd3/qgis-stereonet 

## Install   
Download zip file from github, install into QGIS using plugin manager   

## Usage
 1- Select a layer that has structural info in QGIS   
 2- Select the points you want to plot with one of the Select Tools (**NOT** the Identification Tool)   
 3- You can use the WAXI QFIELD Plugin (https://github.com/swaxi/WAXI_QF) to control display behaviour   
 4- Click on WAXI Stereonet icon    ![plugin_icon](icon.png)  
    
- Planar structures can be displayed as poles or great circles   
- Linear structures are displayed as poles or rose diagrams, but if a planar feature is assocated with a linear feature, that planar feature will optionally be displayed as a great circle in a stereoplot   

## Interactive Stereonet Selection
After plotting a stereonet, you can select poles directly in the plot window and have those features highlighted in the QGIS map layer.

Note: Selection only works in poles mode (i.e., when showGtCircles is false in your config, which is the default). Great-circle and rose-diagram plots do not support selection.

#### Lasso Selection
Click and drag to draw a freehand polygon around any number of poles.

- In the stereonet window, left-click and drag to draw a lasso shape
- Release the mouse button to complete the selection
- All poles inside the lasso are highlighted with red open circles
- The corresponding features are immediately selected in the QGIS map layer
- Shift select adds points to the selection   

#### Single-Point Selection
Click near any individual pole to select it.

- Left-click close to a pole (without dragging)
- The nearest pole within the click tolerance is highlighted
- The corresponding feature is selected in the QGIS map layer
#### Clearing the Selection
Press Escape in the stereonet window to clear all selected poles and remove the selection from the QGIS map layer.

#### Tips
- The lasso and click selection replace the current selection each time — they do not add to an existing selection.
- If a pole is plotted but the click doesn't register, try clicking a little closer to the centre of the point — the tolerance is tuned to the stereonet's projection coordinate space.
- The stereonet window must remain open and in focus for keyboard shortcuts (e.g., Escape) to work.


## Field Names

The following field names are currently recognised, you can go in the file _ _init.py__ from around line 159 and add your own if you like:

- Strike field names = ['Strike_RHR', 'Strike', 'strike']
- Dip Direction field names = ['Dip_Direction', 'Dip_Dir', 'DipDirection', 'dip_direction']
- Dip field names = ['Dip', 'dip']
- Azimuth field names = ['Azimuth', 'azimuth', 'Bearing', 'bearing']
- Plunge field names = ['Plunge', 'plunge']
- Strike of plane for lineations field names = ['Strike_ref', 'Strike_Ref', 'strike_ref']
- Dip of plane for lineations field names = ['Dip_ref', 'Dip_Ref', 'dip_ref']
- Kinematics field names = ['Kinematics', 'kinematics']
- Pitch field names = ['Pitch_RHR', 'Pitch_rhr', 'Pitch_Rhr', 'Pitch', 'pitch_rhr', 'RHR_pitch', 'rhr_pitch', 'pitch']
   
## Roadmap:

-	Plot the hanging wall relative motion vector with an arrow on pole symbol if kinematic information is available.   Broken!
