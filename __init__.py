# -*- coding: utf-8 -*-
"""
/***************************************************************************
 Stereonet
                                 A QGIS plugin
 Displays a geologic stereonet of selected data
                             -------------------
        begin               : 2016-11-29
        copyright           : (C) 2016 by Daniel Childs
        email               : daniel@childsgeo.com
        git sha             : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                          *
 *   This program is free software; you can redistribute it and/or modify   *
 *   it under the terms of the GNU General Public License as published by   *
 *   the Free Software Foundation; either version 2 of the License, or      *
 *   (at your option) any later version.                                    *
 *                                                                          *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""

from qgis.PyQt.QtGui import *
from qgis.PyQt.QtCore import *
from qgis.PyQt.QtWidgets import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from .mplstereonet import *
from qgis.core import *
from qgis.gui import *
import os
from qgis.core import QgsProject
from math import asin,sin,degrees,radians,cos,tan,atan
import json

def classFactory(iface):
    return Stereonet(iface)

class Stereonet:
    def __init__(self, iface):
        self.iface = iface

    def initGui(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self.contourAction = QAction(QIcon(str(dir_path)+"/icon.png"), u'Stereonet', self.iface.mainWindow())
        self.contourAction.triggered.connect(self.contourPlot)
        self.iface.addToolBarIcon(self.contourAction)

    def unload(self):
        self.iface.removeToolBarIcon(self.contourAction)
        del self.contourAction
    
    def waxi_tangent_lineation_plot(self,ax,strikes, dips,kinematics,rhr,azs):
        """Makes a tangent lineation plot for normal faults with the given strikes,
        dips, and rakes."""
        sos=['Sinistral-slip','Dextral-slip','Normal-slip','Reverse-slip']

        for j,sos_it in enumerate(sos):

            dp=list()
            sp=list()
            rh=list()
            az=list()
            for i in range(len(strikes)):
                if(kinematics[i]==sos_it):

                    dp.append(dips[i])
                    sp.append(strikes[i])
                    rh.append(rhr[i])
                    if(sos_it in ['Normal-slip','Reverse-slip']):
                        if(azs[i]>90):
                            azimuth=180-azs[i]
                        else:
                            azimuth=azs[i]
                        if(rhr[i]=='YES'):
                            az.append(360-azimuth)
                        else:
                            az.append(180-azimuth)
                    else:
                        if(azs[i]>180 ):
                            azimuth=90-azs[i]
                        else:
                            azimuth=azs[i]
                        if(rhr[i]=='YES'):
                            azimuth=180-azimuth
                        else:
                            azimuth=360-azimuth
                        if(rhr[i]=='YES' and azs[i]>90):
                            azimuth=azs[i]
                        az.append(azimuth)                        


            # Calculate the position of the rake of the lineations, but don't plot yet
            rake_x, rake_y = mplstereonet.rake(sp, dp, az)
            
            # Calculate the direction the arrows should point
            # These are all normal faults, so the arrows point away from the center
            # Because we're plotting at the pole location, however, we need to flip this
            # from what we plotted with the "ball of string" plot.
            mag = np.hypot(rake_x, rake_y)
            u, v = -rake_x / mag, -rake_y / mag

            # Calculate the position of the poles
            pole_x, pole_y = mplstereonet.pole(sp, dp)
           
            # Plot the arrows centered on the pole locations...
            if(sos_it=='Sinistral-slip'):
                arrows = ax.quiver(pole_x, pole_y, u, v,  width=1, headwidth=4, units='dots', color='r',
                                pivot='tail')
            elif(sos_it=='Dextral-slip'):
                arrows = ax.quiver(pole_x, pole_y, -u, -v,  width=1, headwidth=4, units='dots', color='g',
                                pivot='tail')
            elif(sos_it=='Normal-slip'):
                arrows = ax.quiver(pole_x, pole_y, u,-v,  width=1, headwidth=4, units='dots', color='b',
                                pivot='tail')
            elif(sos_it=='Reverse-slip'):
                arrows = ax.quiver(pole_x, pole_y, -u,v,  width=1, headwidth=4, units='dots', color='m',
                                pivot='tail')


            #return arrows
    


    def rose_diagram(self,strikes,title):
        #modified from: http://geologyandpython.com/structural_geology.html

        bin_edges = np.arange(-5, 366, 10)
        number_of_strikes, bin_edges = np.histogram(strikes, bin_edges)
        number_of_strikes[0] += number_of_strikes[-1]
        half = np.sum(np.split(number_of_strikes[:-1], 2), 0)
        two_halves = np.concatenate([half, half])
        fig = plt.figure(figsize=(8,8))

        """ax = fig.add_subplot(121, projection='stereonet')

        ax.pole(strikes, dips, c='k', label='Pole of the Planes')
        ax.density_contourf(strikes, dips, measurement='poles', cmap='Reds')
        ax.set_title('Density coutour of the Poles', y=1.10, fontsize=15)
        ax.grid()"""

        ax = fig.add_subplot(111, projection='polar')

        ax.bar(np.deg2rad(np.arange(0, 360, 10)), two_halves, 
            width=np.deg2rad(10),  color='.8', edgecolor='k')
        ax.set_theta_zero_location('N')
        ax.set_theta_direction(-1)
        ax.set_thetagrids(np.arange(0, 360, 10), labels=np.arange(0, 360, 10))
        ax.set_rgrids(np.arange(1, two_halves.max() + 1, 2), angle=0, weight= 'black')
        ax.set_title(title)

        fig.tight_layout()
        plt.show()
        
    def contourPlot(self):
        snames=['Strike_RHR', 'Strike', 'strike']
        ddname='Dip_Dir' 
        dname='Dip'
        aname='Azimuth'
        pname='Plunge'
        srefname='Strike_ref'
        drefname='Dip_ref'
        kname='Kinematics'
        pname='Plunge'
        prhrname='Pitch_RHR'

        strikes = list()
        dips = list()
        strikesref = list()
        dipsref = list()
        plunges=list()
        kinematics=list()
        rakes_strikes=list()
        rakes_dips=list()
        rakes_pstrikes=list()
        rakes_pdips=list()
        roseAzimuth = list()
        rhr=list()
        azs=list()


        project = QgsProject.instance()
        proj_file_path=project.fileName()
        head_tail = os.path.split(proj_file_path)
        WAXI_project_path = os.path.abspath(QgsProject.instance().fileName())
        stereoConfigPath = os.path.join(os.path.dirname(WAXI_project_path), "99_COMMAND_FILES_PLUGIN/stereonet.json")

        #stereoConfigPath = head_tail[0]+"/0. FIELD DATA/0. CURRENT MISSION/0. STOPS-SAMPLING-PHOTOGRAPHS-COMMENTS/stereonet.json"
        
        stereoConfig={'showGtCircles':False,'showContours':True,'showKinematics':True,'linPlanes':True,'roseDiagram':False}
        

        if(os.path.exists(stereoConfigPath)):
            with open(stereoConfigPath,"r") as json_file:
                stereoConfig=json.load(json_file)
        
        self.iface.layerTreeView().selectedLayers()

        layers = list(QgsProject.instance().mapLayers().values())
        layers=self.iface.layerTreeView().selectedLayers()

        for layer in layers:
            if layer.type() == QgsMapLayer.VectorLayer:

                iter = layer.selectedFeatures()
                for sn in snames:
                    strikeExists = layer.fields().lookupField(sn)
                    if strikeExists != -1:
                        sname=sn
                        break
                ddrExists = layer.fields().lookupField(ddname)
                dipExists = layer.fields().lookupField(dname)
                azimuthExists = layer.fields().lookupField(aname)
                plungeExists = layer.fields().lookupField(pname)
                srefExists = layer.fields().lookupField(srefname)
                drefExists = layer.fields().lookupField(drefname)

                for feature in iter:
                  
                    if ddrExists != -1 and dipExists != -1:
                        strikes.append(feature[ddname]+90)
                        dips.append(feature[dname])

                    elif strikeExists != -1 and dipExists != -1:
                        strikes.append(feature[sname])
                        dips.append(feature[dname])

                    elif azimuthExists != -1 and plungeExists != -1:
                        strikes.append(feature[aname]+90)
                        dips.append(90-feature[pname])

                    if srefExists != -1 and drefExists != -1:
                        if(not feature[srefname] is None and not feature[drefname] is None):
                            strikesref.append(feature[srefname])
                            dipsref.append(feature[drefname])
                    
                    if plungeExists != -1 and drefExists != -1:
                        if(feature[pname]  and  feature[drefname]  and  feature[kname] ):
                            rakes_strikes.append(feature[srefname])
                            rakes_dips.append(feature[drefname])
                            kinematics.append(feature[kname])
                            rhr.append(feature[prhrname])
                            azs.append(feature[aname])

                    if azimuthExists != -1 and stereoConfig['roseDiagram']:
                        roseAzimuth.append(feature[aname])
 


            else:
                continue

        strikesref = [i for i in strikesref if i != None]
        dipsref = [i for i in dipsref if i != None]
        strikes = [i for i in strikes if i != None]
        dips = [i for i in dips if i != None]
        plunges = [i for i in plunges if i != None]

        if(len(roseAzimuth) != 0 and stereoConfig['roseDiagram']):
            self.rose_diagram(roseAzimuth,layer.name()+" [# "+str(len(iter))+"]")
        elif (len(strikes) != 0):
            fig, ax = mplstereonet.subplots()
            ax.set_azimuth_ticks([0,30,60,90,120,150,180,210,240,270,300,330])
            ax.set_azimuth_ticklabels(['0\u00b0','30\u00b0','60\u00b0','90\u00b0','120\u00b0','150\u00b0','180\u00b0','210\u00b0','240\u00b0','270\u00b0','300\u00b0','330\u00b0'])
            ax.grid(kind='equal_area_stereonet')
            if(stereoConfig['showContours']):
                ax.density_contour(strikes, dips, measurement='poles',cmap=cm.coolwarm,method='exponential_kamb',sigma=1.5,linewidths =0.5)
            if(stereoConfig['showGtCircles'] and strikeExists != -1):
                ax.plane(strikes, dips, 'k',linewidth=1)
            else:
                ax.pole(strikes, dips, 'k.', markersize=5)
                if(srefExists != -1 and drefExists != -1 and stereoConfig['linPlanes']):
                    ax.plane(strikesref,dipsref,'k',linewidth=1)
                if plungeExists != -1 and drefExists != -1 and stereoConfig['showKinematics']:
                    self.waxi_tangent_lineation_plot(ax,rakes_strikes, rakes_dips,kinematics,rhr,azs)

            ax.set_title(layer.name()+" [# "+str(len(iter))+"]",pad=24)
            plt.show()

        else:
            self.iface.messageBar().pushMessage("No data selected, or no structural data found: first select a layer with structural info, then select the points that you wish to plot", level=Qgis.Warning, duration=5)