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
   
## Roadmap:

-	Plot the hanging wall relative motion vector with an arrow on pole symbol if kinematic information is available.   Broken!
