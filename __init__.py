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

from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from .mplstereonet import *
from qgis.core import *
from qgis.gui import *
import os
from qgis.core import QgsProject
from math import asin,sin,degrees,radians,cos
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

    def plungedip2rake(self,plunge,dip,kinematics):
        if(plunge==dip):
            rake=0
        elif(dip==90.0):
            rake=90-plunge
        else:
            rake=90-degrees(asin(sin(radians(plunge))/sin(radians(dip))))
        
        return(rake)

    def waxi_fault_and_striae_plot(self,ax, pstrikes,pdips,strikes, dips, rakes,kinematics):
        sos=['Sinistral-slip','Dextral-Slip','Normal-slip','Reverse-slip']
        """Makes a fault-and-striae plot (a.k.a. "Ball of String") for normal faults
        with the given strikes, dips, and rakes."""
        # Plot the planes
        #print(kinematics)
        #print(strikes)
        #print(dips)
        #print(rakes)

        # loop through each kinematic sense to plot correct orientations of arrows
        for sos_it in sos:
            s=list()
            d=list()
            r=list()            
            dp=list()
            sp=list()
            for i in range(len(strikes)):
                if(kinematics[i]==sos_it):
                    s.append(pstrikes[i])
                    d.append(pdips[i])
                    r.append(rakes[i])
                    dp.append(dips[i])
                    sp.append(strikes[i])
                    #print(i,strikes[i],dips[i],rakes[i])
            # Calculate the position of the rake of the lineations, but don't plot yet
            x,y = mplstereonet.rake(s, d,r)
            x1,y1 = mplstereonet.pole(sp, dp)

            # Calculate the direction the arrows should point
            # These are all normal faults, so the arrows point away from the center
            # For thrusts, it would just be u, v = -x/mag, -y/mag
            mag = np.hypot(x, y)

            theta = np.deg2rad(0)

            rot = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])

            v1 = np.array([x, y])

            v2 = np.matmul(rot, v1)
            u, v = v2[0] / mag, v2[1] / mag

            # Plot the arrows at the rake locations...
            #print(sos_it,s,d,r)
            if(sos_it == sos[0] and len(x)>0):
                #print("rot")
                #print(v)
                #print(v2)

                arrows = ax.quiver(x1, y1, v, u, width=1, headwidth=4, units='dots',color='r')
                #arrows = ax.quiver(x1, y1, -v, -u, width=1, headwidth=4, units='dots',color='b')
                #arrows = ax.quiver(x1, y1, u, -v, width=1, headwidth=4, units='dots',color='g')
                #arrows = ax.quiver(x1, y1, -u, v, width=1, headwidth=4, units='dots',color='k')
            elif(sos_it == sos[1]and len(x)>0):
                arrows = ax.quiver(x1, y1, -v, -u, width=1, headwidth=4, units='dots',color='g')
            elif(sos_it == sos[2]and len(x)>0):
                arrows = ax.quiver(x1, y1, -u, v, width=1, headwidth=4, units='dots',color='b')
            elif(sos_it == sos[3]and len(x)>0):
                #arrows = ax.quiver(x1, y1, v, u, width=1, headwidth=4, units='dots',color='r')
                #arrows = ax.quiver(x1, y1, -v, -u, width=1, headwidth=4, units='dots',color='b')
                arrows = ax.quiver(x1, y1, u, -v, width=1, headwidth=4, units='dots',color='k')
                #arrows = ax.quiver(x1, y1, -u, v, width=1, headwidth=4, units='dots',color='k')


    def contourPlot(self):
        sname='Strike_RHR'
        ddname='Dip_Dir'
        dname='Dip'
        aname='Azimuth'
        pname='Plunge'
        srefname='Strike_ref'
        drefname='Dip_ref'
        kname='Kinematics'
        pname='Plunge'

        strikes = list()
        dips = list()
        strikesref = list()
        dipsref = list()
        plunges=list()
        kinematics=list()
        rakes=list()
        rakes_strikes=list()
        rakes_dips=list()
        rakes_pstrikes=list()
        rakes_pdips=list()
        project = QgsProject.instance()
        proj_file_path=project.fileName()
        head_tail = os.path.split(proj_file_path)

        stereoConfigPath = head_tail[0]+"/0. FIELD DATA/0. CURRENT MISSION/0. STOPS-SAMPLING-PHOTOGRAPHS-COMMENTS/stereonet.json"
        
        stereoConfig={'showGtCircles':True,'showContours':True,'showKinematics':True}
        

        if(os.path.exists(stereoConfigPath)):
            with open(stereoConfigPath,"r") as json_file:
                stereoConfig=json.load(json_file)
        
        self.iface.layerTreeView().selectedLayers()

        layers = list(QgsProject.instance().mapLayers().values())
        layers=self.iface.layerTreeView().selectedLayers()

        for layer in layers:
            if layer.type() == QgsMapLayer.VectorLayer:
                iter = layer.selectedFeatures()
                strikeExists = layer.fields().lookupField(sname)
                ddrExists = layer.fields().lookupField(ddname)
                dipExists = layer.fields().lookupField(dname)
                azimuthExists = layer.fields().lookupField(aname)
                plungeExists = layer.fields().lookupField(pname)
                srefExists = layer.fields().lookupField(srefname)
                drefExists = layer.fields().lookupField(drefname)

                for feature in iter:
                    if strikeExists != -1 and dipExists != -1:
                        strikes.append(feature[sname])
                        dips.append(feature[dname])
                  
                    elif ddrExists != -1 and dipExists != -1:
                        strikes.append(feature[ddname]+90)
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
                            rakes.append(self.plungedip2rake(feature[pname],feature[drefname],feature[kname]))
                            rakes_strikes.append(feature[srefname])
                            rakes_dips.append(feature[drefname])
                            kinematics.append(feature[kname])
                            rakes_pstrikes.append(feature[aname]+90)
                            rakes_pdips.append(90-feature[pname])


            else:
                continue
        strikesref = [i for i in strikesref if i != None]
        dipsref = [i for i in dipsref if i != None]
        strikes = [i for i in strikes if i != None]
        dips = [i for i in dips if i != None]
        plunges = [i for i in plunges if i != None]
        #print(strikes)

        if (len(strikes) != 0):
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
                if(srefExists != -1 and drefExists != -1 ):
                    ax.plane(strikesref,dipsref,'k',linewidth=1)
                if plungeExists != -1 and drefExists != -1 and stereoConfig['showKinematics']:
                    self.waxi_fault_and_striae_plot(ax, rakes_pstrikes, rakes_pdips,rakes_strikes, rakes_dips, rakes,kinematics)

            ax.set_title(layer.name()+" [# "+str(len(iter))+"]")
            plt.show()
        else:
            self.iface.messageBar().pushMessage("No data selected, or no structural data found: first select a layer with structural info, then select the points that you wish to plot", level=Qgis.Warning, duration=5)