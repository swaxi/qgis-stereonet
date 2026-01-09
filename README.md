# qgis-stereonet
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
