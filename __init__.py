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

    def contourPlot(self):
        sname='Strike_RHR'
        ddname='Dip_Dir'
        dname='Dip'
        aname='Azimuth'
        pname='Plunge'
        strikes = list()
        dips = list()
 #       layers = self.iface.legendInterface().layers()

        
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

            else:
                continue

        if len(strikes) != 0:
            fig, ax = mplstereonet.subplots()
            ax.set_azimuth_ticks([0,30,60,90,120,150,180,210,240,270,300,330])
            ax.set_azimuth_ticklabels(['0\u00b0','30\u00b0','60\u00b0','90\u00b0','120\u00b0','150\u00b0','180\u00b0','210\u00b0','240\u00b0','270\u00b0','300\u00b0','330\u00b0'])
            ax.grid(kind='equal_area_stereonet')
            ax.density_contour(strikes, dips, measurement='poles',cmap=cm.coolwarm,method='exponential_kamb',sigma=1.5,linewidths =0.5)
            ax.pole(strikes, dips, 'k.', markersize=7)
            ax.set_title(layer.name()+" [# "+str(len(iter))+"]")
            plt.show()
        else:
            self.iface.messageBar().pushMessage("No data selected, or no structural data found: first select a layer with structural info, then select the points that you wish to plot", level=Qgis.Warning, duration=5)