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
from matplotlib.path import Path
from collections import defaultdict
from .mplstereonet import *
from qgis.core import *
from qgis.gui import *
import os
from qgis.core import QgsProject
from math import asin,sin,degrees,radians,cos,tan,atan
import json

class _BoundedLassoSelector:
    """Lasso selector that works entirely in axes coordinates (0-1 space).

    This avoids the numerically unstable inverse Lambert projection transform.
    The stereonet boundary is always the circle centred at (0.5, 0.5) with
    radius 0.5 in axes coordinates.  When the cursor leaves that circle the
    path is clipped to the boundary; arc vertices are interpolated so the line
    visually hugs the edge.  The onselect callback receives a list of
    (ax_x, ax_y) vertices in axes coordinates.
    """

    _R = 0.5  # stereonet radius in axes coords

    def __init__(self, ax, onselect):
        self.ax = ax
        self.onselect = onselect
        self._active = False
        self._verts = []          # list of (ax_x, ax_y)
        self._last_clipped = False
        self._last_angle = 0.0
        # Line drawn in axes-coordinate space so no data transform is needed
        self._line, = ax.plot([], [], color='black', linewidth=0.8,
                              transform=ax.transAxes)
        self._line.set_visible(False)
        canvas = ax.figure.canvas
        self._cids = [
            canvas.mpl_connect('button_press_event', self._on_press),
            canvas.mpl_connect('button_release_event', self._on_release),
            canvas.mpl_connect('motion_notify_event', self._on_move),
        ]

    def _disp_to_axes(self, x, y):
        """Convert a single display-coord point to axes coords."""
        return self.ax.transAxes.inverted().transform([[x, y]])[0]

    def _clip(self, x, y):
        """Convert display coords to axes coords, clipped to the stereonet circle.
        Returns (ax_x, ax_y, is_clipped, angle_from_centre)."""
        ax_x, ax_y = self._disp_to_axes(x, y)
        dx, dy = ax_x - 0.5, ax_y - 0.5
        angle = np.arctan2(dy, dx)
        if np.hypot(dx, dy) <= self._R:
            return ax_x, ax_y, False, angle
        return (0.5 + self._R * np.cos(angle),
                0.5 + self._R * np.sin(angle),
                True, angle)

    def _arc_verts(self, from_angle, to_angle):
        """Axes-coord vertices on the boundary arc, start-exclusive, end-inclusive."""
        diff = (to_angle - from_angle + np.pi) % (2 * np.pi) - np.pi
        n = max(2, int(abs(diff) * 20))
        angles = np.linspace(from_angle, from_angle + diff, n + 1)[1:]
        return [(0.5 + self._R * np.cos(a), 0.5 + self._R * np.sin(a))
                for a in angles]

    def _update_line(self):
        if self._verts:
            xs, ys = zip(*self._verts)
        else:
            xs, ys = [], []
        self._line.set_data(xs, ys)

    def _on_press(self, event):
        if event.button != 1:
            return
        bbox = self.ax.get_window_extent()
        if not (bbox.x0 <= event.x <= bbox.x1 and bbox.y0 <= event.y <= bbox.y1):
            return
        self._active = True
        self._verts = []
        ax_x, ax_y, clipped, angle = self._clip(event.x, event.y)
        self._verts.append((ax_x, ax_y))
        self._last_clipped = clipped
        self._last_angle = angle
        self._line.set_visible(True)
        self._update_line()
        self.ax.figure.canvas.draw_idle()

    def _on_move(self, event):
        if not self._active:
            return
        ax_x, ax_y, clipped, angle = self._clip(event.x, event.y)
        if self._last_clipped and clipped:
            self._verts.extend(self._arc_verts(self._last_angle, angle))
            self._last_angle = angle
        else:
            self._verts.append((ax_x, ax_y))
            self._last_clipped = clipped
            if clipped:
                self._last_angle = angle
        self._update_line()
        self.ax.figure.canvas.draw_idle()

    def _on_release(self, event):
        if not self._active or event.button != 1:
            return
        self._active = False
        self._line.set_visible(False)
        self.ax.figure.canvas.draw_idle()
        self.onselect(self._verts)

    def disconnect(self):
        for cid in self._cids:
            self.ax.figure.canvas.mpl_disconnect(cid)


class StereonetSettingsDialog(QDialog):
    """Settings dialog for controlling stereonet plot style.

    If config_path points to an existing stereonet.json the dialog reads from
    and writes back to that file.  Otherwise values are persisted in QSettings.
    """

    _QSETTINGS_ORG = 'qgis-stereonet'
    _QSETTINGS_APP = 'stereonet'
    _DEFAULTS = {
        'showGtCircles': False, 'showContours': True,
        'showKinematics': True, 'linPlanes': True, 'roseDiagram': False,
        'fitGirdle': False, 'dataType': 'Planes Only',
    }

    def __init__(self, parent=None, config_path=None, detected_data_type=None):
        super().__init__(parent)
        self._config_path = config_path
        self.setWindowTitle('Stereographic Projection Settings')
        self.setModal(True)

        cfg = self._load()

        outer = QVBoxLayout()
        outer.addWidget(QLabel(
            'Select Features to Plot '
            '(Lower Hemisphere, Equal-Area Stereonet Projection):'))

        row = QHBoxLayout()
        self.gtCircles_cb  = QCheckBox('Great Circles')
        self.contours_cb   = QCheckBox('Contours')
        self.linPlanes_cb  = QCheckBox('Lineation-bearing Planes')
        self.rose_cb       = QCheckBox('Rose Diagram')
        self.kinematics_cb = QCheckBox('Kinematics')
        self.fitGirdle_cb  = QCheckBox('Best Fit Girdle')

        self.gtCircles_cb.setChecked( cfg['showGtCircles'])
        self.contours_cb.setChecked(  cfg['showContours'])
        self.linPlanes_cb.setChecked( cfg['linPlanes'])
        self.rose_cb.setChecked(      cfg['roseDiagram'])
        self.kinematics_cb.setChecked(cfg['showKinematics'])
        self.fitGirdle_cb.setChecked( cfg['fitGirdle'])

        for cb in [self.gtCircles_cb, self.contours_cb, self.linPlanes_cb,
                   self.rose_cb, self.kinematics_cb, self.fitGirdle_cb]:
            row.addWidget(cb)

        update_btn = QPushButton('Update Settings')
        update_btn.clicked.connect(self._save_and_close)
        row.addStretch()
        row.addWidget(update_btn)

        outer.addLayout(row)

        dt_row = QHBoxLayout()
        dt_row.addWidget(QLabel('Data to plot:'))
        self.dataType_cb = QComboBox()
        self.dataType_cb.addItems([
            'Planes Only',
            'Lineations Only',
            'Lineations with Bearing Planes'
        ])

        initial_data_type = detected_data_type or cfg.get('dataType', 'Planes Only')
        if initial_data_type == 'Lineations with Planes':
            initial_data_type = 'Lineations with Bearing Planes'
        self.dataType_cb.setCurrentText(initial_data_type)

        self.dataType_cb.currentTextChanged.connect(self._on_data_type_changed)
        self._on_data_type_changed(self.dataType_cb.currentText())
        dt_row.addWidget(self.dataType_cb)
        dt_row.addStretch()
        outer.addLayout(dt_row)

        self.setLayout(outer)

    def _on_data_type_changed(self, text):
        """Synchronise the bearing-plane checkbox with the selected data mode."""
        if text == 'Lineations with Bearing Planes':
            self.linPlanes_cb.setChecked(True)
            self.linPlanes_cb.setEnabled(True)
        else:
            self.linPlanes_cb.setChecked(False)
            self.linPlanes_cb.setEnabled(False)

    def _load(self):
        """Load config from JSON file if present, else QSettings, else defaults."""
        if self._config_path and os.path.exists(self._config_path):
            with open(self._config_path, 'r') as f:
                cfg = json.load(f)
            return {k: cfg.get(k, v) for k, v in self._DEFAULTS.items()}
        s = QSettings(self._QSETTINGS_ORG, self._QSETTINGS_APP)
        if s.contains('showContours'):
            result = {}
            for k, v in self._DEFAULTS.items():
                result[k] = s.value(k, v, type=bool) if isinstance(v, bool) else s.value(k, v)
            return result
        return dict(self._DEFAULTS)

    def _save_and_close(self):
        if self.dataType_cb.currentText() == 'Lineations with Bearing Planes':
            self.linPlanes_cb.setChecked(True)
        else:
            self.linPlanes_cb.setChecked(False)

        cfg = {
            'showGtCircles':  self.gtCircles_cb.isChecked(),
            'showContours':   self.contours_cb.isChecked(),
            'showKinematics': self.kinematics_cb.isChecked(),
            'linPlanes':      self.linPlanes_cb.isChecked(),
            'roseDiagram':    self.rose_cb.isChecked(),
            'fitGirdle':      self.fitGirdle_cb.isChecked(),
            'dataType':       self.dataType_cb.currentText(),
        }
        if self._config_path:
            os.makedirs(os.path.dirname(self._config_path), exist_ok=True)
            with open(self._config_path, 'w') as f:
                json.dump(cfg, f, indent=4)
        else:
            s = QSettings(self._QSETTINGS_ORG, self._QSETTINGS_APP)
            for k, v in cfg.items():
                s.setValue(k, v)
        self.accept()

    @staticmethod
    def load_qsettings():
        """Return a stereoConfig dict from QSettings, or None if not yet saved."""
        s = QSettings(StereonetSettingsDialog._QSETTINGS_ORG,
                      StereonetSettingsDialog._QSETTINGS_APP)
        if not s.contains('showContours'):
            return None
        result = {}
        for k, v in StereonetSettingsDialog._DEFAULTS.items():
            result[k] = s.value(k, v, type=bool) if isinstance(v, bool) else s.value(k, v)
        return result


def _attr(v):
    """Return a QGIS feature attribute as a Python value, or None if null."""
    if v is None:
        return None
    if hasattr(v, 'isNull') and v.isNull():
        return None
    return v


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

        self.settingsAction = QAction(
            QgsApplication.getThemeIcon('mActionOptions.svg'),
            u'Stereonet Settings', self.iface.mainWindow())
        self.settingsAction.triggered.connect(self.showSettings)
        self.iface.addToolBarIcon(self.settingsAction)

    def unload(self):
        self.iface.removeToolBarIcon(self.contourAction)
        del self.contourAction
        self.iface.removeToolBarIcon(self.settingsAction)
        del self.settingsAction

    def _structural_field_names(self):
        """Return accepted field-name variants for structural data detection."""
        return {
            'strike': ['Strike_RHR', 'Strike', 'strike'],
            'dipdir': ['Dip_Direction', 'Dip_Dir', 'DipDirection',
                       'dip_direction', 'DipDir', 'DIPDIR'],
            'dip': ['Dip', 'dip'],
            'azimuth': ['Azimuth', 'azimuth', 'Bearing', 'bearing',
                        'Trend', 'TREND'],
            'plunge': ['Plunge', 'plunge'],
            'strike_ref': ['Strike_ref', 'Strike_Ref', 'strike_ref'],
            'dip_ref': ['Dip_ref', 'Dip_Ref', 'dip_ref'],
            'dipdir_ref': ['DipDir_ref', 'Dip_Dir_ref', 'DipDirection_ref',
                           'Dip_Direction_ref', 'DipDirection_Ref',
                           'dipdir_ref', 'dip_dir_ref', 'dip_direction_ref'],
            'kinematics': ['Kinematics', 'kinematics'],
            'pitch_rhr': ['Pitch_RHR', 'Pitch_rhr', 'Pitch_Rhr', 'Pitch',
                          'pitch_rhr', 'RHR_pitch', 'rhr_pitch', 'pitch'],
        }

    def _field_exists(self, layer, fieldnames):
        for fieldname in fieldnames:
            field_index = layer.fields().lookupField(fieldname)
            if field_index != -1:
                return True, fieldname
        return False, False

    def _detect_data_type_from_layers(self, layers):
        """Infer the plotting mode from fields available in selected layers."""
        names = self._structural_field_names()
        has_plane = False
        has_line = False
        has_ref_plane = False

        for layer in layers:
            if layer.type() != QgsMapLayer.VectorLayer:
                continue

            strike_ok, _ = self._field_exists(layer, names['strike'])
            dipdir_ok, _ = self._field_exists(layer, names['dipdir'])
            dip_ok, _ = self._field_exists(layer, names['dip'])
            az_ok, _ = self._field_exists(layer, names['azimuth'])
            plunge_ok, _ = self._field_exists(layer, names['plunge'])
            sref_ok, _ = self._field_exists(layer, names['strike_ref'])
            ddref_ok, _ = self._field_exists(layer, names['dipdir_ref'])
            dref_ok, _ = self._field_exists(layer, names['dip_ref'])

            has_plane = has_plane or ((strike_ok or dipdir_ok) and dip_ok)
            has_line = has_line or (az_ok and plunge_ok)
            has_ref_plane = has_ref_plane or ((sref_ok or ddref_ok) and dref_ok)

        # Combined mode applies to linear data carrying either explicit
        # reference-plane fields (e.g. Strike_ref/Dip_ref or DipDir_ref/Dip_ref)
        # or regular planar fields on the same layer (e.g. Folds_PT with
        # Azimuth/Plunge plus Strike_RHR or Dip_Dir/Dip).
        if has_line and (has_ref_plane or has_plane):
            return 'Lineations with Bearing Planes'
        if has_line:
            return 'Lineations Only'
        if has_plane:
            return 'Planes Only'
        return None

    def showSettings(self):
        project_file = QgsProject.instance().fileName()
        config_path = None
        if project_file:
            config_path = os.path.join(
                os.path.dirname(os.path.abspath(project_file)),
                "99_COMMAND_FILES_PLUGIN",
                "stereonet.json")

        layers = self.iface.layerTreeView().selectedLayers()
        detected_data_type = self._detect_data_type_from_layers(layers)

        dlg = StereonetSettingsDialog(
            self.iface.mainWindow(),
            config_path=config_path,
            detected_data_type=detected_data_type)
        dlg.exec()
    
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
        names = self._structural_field_names()
        snames = names['strike']
        ddnames = names['dipdir']
        dnames = names['dip']
        anames = names['azimuth']
        pnames = names['plunge']
        srefnames = names['strike_ref']
        drefnames = names['dip_ref']
        ddrefnames = names['dipdir_ref']
        knames = names['kinematics']
        prhrnames = names['pitch_rhr']
        
        plane_strikes = list()
        plane_dips = list()
        plane_feature_ids = list()
        plane_labels = list()
        linear_plunges = list()
        linear_bearings = list()
        linear_feature_ids = list()
        linear_labels = list()
        strikesref = list()
        dipsref = list()
        plunges = list()
        kinematics = list()
        rakes_strikes = list()
        rakes_dips = list()
        roseAzimuth = list()
        rhr = list()
        azs = list()


        project = QgsProject.instance()
        proj_file_path=project.fileName()
        head_tail = os.path.split(proj_file_path)
        WAXI_project_path = os.path.abspath(QgsProject.instance().fileName())
        stereoConfigPath = os.path.join(os.path.dirname(WAXI_project_path), "99_COMMAND_FILES_PLUGIN/stereonet.json")

        #stereoConfigPath = head_tail[0]+"/0. FIELD DATA/0. CURRENT MISSION/0. STOPS-SAMPLING-PHOTOGRAPHS-COMMENTS/stereonet.json"
        
        stereoConfig = {'showGtCircles': False, 'showContours': True,
                        'showKinematics': True, 'linPlanes': True, 'roseDiagram': False,
                        'fitGirdle': False, 'dataType': 'Planes Only'}

        if os.path.exists(stereoConfigPath):
            with open(stereoConfigPath, "r") as json_file:
                stereoConfig = json.load(json_file)
        else:
            qs = StereonetSettingsDialog.load_qsettings()
            if qs is not None:
                stereoConfig = qs
        
        self.iface.layerTreeView().selectedLayers()

        layers = list(QgsProject.instance().mapLayers().values())
        layers=self.iface.layerTreeView().selectedLayers()

        for layer in layers:
            if layer.type() == QgsMapLayer.VectorLayer:

                iter = layer.selectedFeatures()
                strikeExists, sname = self._field_exists(layer,snames)
                ddrExists, ddname = self._field_exists(layer,ddnames)
                dipExists, dname = self._field_exists(layer,dnames)
                azimuthExists, aname = self._field_exists(layer,anames)
                plungeExists, pname = self._field_exists(layer,pnames)
                srefExists, srefname = self._field_exists(layer,srefnames)
                drefExists, drefname = self._field_exists(layer,drefnames)
                ddrefExists, ddrefname = self._field_exists(layer,ddrefnames)
                kinematicsExists, kname = self._field_exists(layer,knames)
                prhrExists, prhrname = self._field_exists(layer,prhrnames )




                for feature in iter:
                    # Capture plane data (dip direction/dip or strike/dip)
                    if ddrExists and dipExists:
                        val_dd, val_d = _attr(feature[ddname]), _attr(feature[dname])
                        if val_dd is not None and val_d is not None:
                            plane_strikes.append((int(val_dd) - 90) % 360)
                            plane_dips.append(float(val_d))
                            plane_feature_ids.append((layer, feature.id()))
                            plane_labels.append(f"{int(val_d)}/{int(val_dd):03d}")
                    elif strikeExists and dipExists:
                        val_s, val_d = _attr(feature[sname]), _attr(feature[dname])
                        if val_s is not None and val_d is not None:
                            plane_strikes.append(float(val_s))
                            plane_dips.append(float(val_d))
                            plane_feature_ids.append((layer, feature.id()))
                            plane_labels.append(f"{int(val_d)}/{int(val_s):03d}")

                    # Capture linear data (azimuth/plunge) independently
                    if azimuthExists and plungeExists:
                        val_a, val_p = _attr(feature[aname]), _attr(feature[pname])
                        if val_a is not None and val_p is not None:
                            linear_plunges.append(int(val_p))
                            linear_bearings.append(int(val_a))
                            linear_feature_ids.append((layer, feature.id()))
                            linear_labels.append(f"{int(val_p)}/{int(val_a):03d}")

                    if drefExists:
                        vd = _attr(feature[drefname])
                        if vd is not None and srefExists:
                            vs = _attr(feature[srefname])
                            if vs is not None:
                                strikesref.append(float(vs))
                                dipsref.append(float(vd))
                        elif vd is not None and ddrefExists:
                            vdd = _attr(feature[ddrefname])
                            if vdd is not None:
                                strikesref.append((float(vdd) - 90) % 360)
                                dipsref.append(float(vd))

                    if (plungeExists and drefExists and
                            kinematicsExists and azimuthExists):
                        vp = _attr(feature[pname])
                        vdr = _attr(feature[drefname])
                        vk = _attr(feature[kname])
                        if vp is not None and vdr is not None and vk is not None:
                            if srefExists:
                                rake_strike = _attr(feature[srefname])
                            elif ddrefExists and _attr(feature[ddrefname]) is not None:
                                rake_strike = (float(_attr(feature[ddrefname])) - 90) % 360
                            else:
                                rake_strike = None
                            if rake_strike is not None:
                                rakes_strikes.append(rake_strike)
                                rakes_dips.append(vdr)
                                kinematics.append(vk)
                                rhr.append(_attr(feature[prhrname]) if prhrExists else None)
                                azs.append(_attr(feature[aname]))

                    if azimuthExists and stereoConfig.get('roseDiagram', False):
                        va = _attr(feature[aname])
                        if va is not None:
                            roseAzimuth.append(va)
 


            else:
                continue

        strikesref = [i for i in strikesref if i is not None]
        dipsref = [i for i in dipsref if i is not None]

        # Determine effective data type from the selected layer fields, then
        # fall back to saved settings if detection is inconclusive.
        detected_data_type = self._detect_data_type_from_layers(layers)
        effective_data_type = detected_data_type or stereoConfig.get('dataType', 'Planes Only')
        if effective_data_type == 'Lineations with Planes':
            effective_data_type = 'Lineations with Bearing Planes'

        show_planes = effective_data_type == 'Planes Only'
        show_linears = effective_data_type in (
            'Lineations Only',
            'Lineations with Bearing Planes')
        show_bearing_planes = effective_data_type == 'Lineations with Bearing Planes'
        # In combined lineation/reference-plane mode, the bearing planes are
        # part of the selected data type itself.  Do not rely on a stale saved
        # linPlanes value from QSettings/stereonet.json; otherwise the checkbox
        # can appear ON in the dialog while the first plot still suppresses the
        # reference-plane great circles until the user toggles and re-saves it.
        effective_lin_planes = (
            True if show_bearing_planes
            else stereoConfig.get('linPlanes', True)
        )
        has_planes = len(plane_strikes) > 0
        has_linears = len(linear_plunges) > 0
        has_bearing_planes = len(strikesref) > 0 and len(dipsref) > 0

        if len(roseAzimuth) != 0 and stereoConfig.get('roseDiagram', False):
            self.rose_diagram(roseAzimuth, layer.name() + " [# " + str(len(iter)) + "]")
        elif ((show_planes and has_planes) or
              (show_linears and has_linears) or
              (show_bearing_planes and has_bearing_planes)):
            fig, ax = mplstereonet.subplots()
            ax.set_azimuth_ticks([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330])
            ax.set_azimuth_ticklabels(['0\u00b0', '30\u00b0', '60\u00b0', '90\u00b0',
                                        '120\u00b0', '150\u00b0', '180\u00b0', '210\u00b0',
                                        '240\u00b0', '270\u00b0', '300\u00b0', '330\u00b0'])
            ax.grid(kind='equal_area_stereonet')

            pole_lines = None
            lin_lines = None

            if show_planes and has_planes:
                if stereoConfig.get('showContours', True):
                    ax.density_contour(plane_strikes, plane_dips, measurement='poles',
                                       cmap=cm.coolwarm, method='exponential_kamb',
                                       sigma=1.5, linewidths=0.5)
                if stereoConfig.get('showGtCircles', False):
                    ax.plane(plane_strikes, plane_dips, 'k', linewidth=1)
                else:
                    pole_lines = ax.pole(plane_strikes, plane_dips, 'k.', markersize=5)

                if stereoConfig.get('fitGirdle', False) and len(plane_strikes) >= 3:
                    gs, gd = mplstereonet.fit_girdle(plane_strikes, plane_dips,
                                                     measurement='poles')
                    ax.plane(gs, gd, 'b-', linewidth=1.5)
                    ax.pole(gs, gd, 'ro', markersize=8)
                    plunge, bearing = mplstereonet.pole2plunge_bearing(gs, gd)
                    _, _, evals = mplstereonet.eigenvectors(
                        plane_strikes, plane_dips, measurement='poles')
                    e1, e2, e3 = evals[0], evals[1], evals[2]
                    pb = (f'Pole to best fit girdle\n'
                          f'  Plunge/Bearing: '
                          f'{int(round(plunge[0]))}/{int(round(bearing[0])):03d}')
                    if e1 > e2 > e3 > 1e-10:
                        K = np.log(e1 / e2) / np.log(e2 / e3)
                        C = np.log(e1 / e3)
                        kshape = ('girdle' if K < 0.9 else
                                  'cluster' if K > 1.1 else 'transitional')
                        info = (
                            f'{pb}\n'
                            f'  K = {K:.2f} ({kshape}),  C = {C:.2f}\n'
                            f'{"─" * 32}\n'
                            f'K (shape):\n'
                            f'  <1 = girdle  ≈1 = transitional  >1 = cluster\n'
                            f'C (strength):\n'
                            f'  0 = random → larger = stronger fabric'
                        )
                    else:
                        info = pb
                    ax.text(.9, -0.1, info,
                            transform=ax.transAxes, ha='left', va='bottom',
                            fontsize=7, clip_on=False,
                            bbox=dict(boxstyle='round,pad=0.4', fc='lightblue',
                                      ec='steelblue', alpha=0.85))

            if show_linears and has_linears:
                if stereoConfig.get('showContours', True):
                    ax.density_contour(linear_plunges, linear_bearings,
                                       measurement='lines', cmap=cm.coolwarm,
                                       method='exponential_kamb', sigma=1.5,
                                       linewidths=0.5)
                lin_lines = ax.line(linear_plunges, linear_bearings,
                                    'k.', markersize=5)

            if show_bearing_planes and effective_lin_planes:
                # Prefer explicit reference-plane fields on linear layers. If
                # absent, use the regular planar fields carried by the same
                # selected layer, as in Folds_PT datasets with Dip/DipDir or
                # Strike plus Plunge/Trend.
                if has_bearing_planes:
                    ax.plane(strikesref, dipsref, 'k', linewidth=1)
                elif has_planes:
                    ax.plane(plane_strikes, plane_dips, 'k', linewidth=1)

                if stereoConfig.get('showKinematics', True) and rakes_strikes:
                    self.waxi_tangent_lineation_plot(ax, rakes_strikes, rakes_dips,
                                                     kinematics, rhr, azs)

            # Resolve which plotted points to use for interactive selection
            pts = None
            sel_fids = None
            sel_labels = None
            if lin_lines and linear_feature_ids:
                lon_data = np.array(lin_lines[0].get_xdata())
                lat_data = np.array(lin_lines[0].get_ydata())
                valid_pts = ~(np.isnan(lon_data) | np.isnan(lat_data))
                if valid_pts.any():
                    pts = np.column_stack([lon_data[valid_pts], lat_data[valid_pts]])
                    sel_fids   = [linear_feature_ids[i] for i, ok in enumerate(valid_pts) if ok]
                    sel_labels = [linear_labels[i]      for i, ok in enumerate(valid_pts) if ok]
            elif pole_lines and plane_feature_ids:
                lon_data = np.array(pole_lines[0].get_xdata())
                lat_data = np.array(pole_lines[0].get_ydata())
                valid_pts = ~(np.isnan(lon_data) | np.isnan(lat_data))
                if valid_pts.any():
                    pts = np.column_stack([lon_data[valid_pts], lat_data[valid_pts]])
                    sel_fids   = [plane_feature_ids[i] for i, ok in enumerate(valid_pts) if ok]
                    sel_labels = [plane_labels[i]      for i, ok in enumerate(valid_pts) if ok]

            if pts is not None and len(pts) > 0:
                sel_plot, = ax.plot([], [], 'ro', markersize=10, zorder=5,
                                    fillstyle='none', markeredgewidth=2)

                annot = ax.annotate(
                    '', xy=(0, 0), xycoords='data',
                    xytext=(8, 8), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.3', fc='lightyellow',
                              ec='gray', alpha=0.9),
                    fontsize=9, zorder=10)
                annot.set_visible(False)
                _hover_idx = [-1]

                def _on_hover(event):
                    if len(pts) == 0:
                        return
                    if event.inaxes != ax:
                        if _hover_idx[0] >= 0:
                            _hover_idx[0] = -1
                            annot.set_visible(False)
                            fig.canvas.draw_idle()
                        return
                    pts_axes = ax.transAxes.inverted().transform(
                        ax.transData.transform(pts))
                    cur_axes = ax.transAxes.inverted().transform(
                        [[event.x, event.y]])[0]
                    dists = np.hypot(pts_axes[:, 0] - cur_axes[0],
                                     pts_axes[:, 1] - cur_axes[1])
                    min_i = int(np.argmin(dists))
                    if dists[min_i] < 0.04:
                        if min_i != _hover_idx[0]:
                            _hover_idx[0] = min_i
                            annot.xy = (pts[min_i, 0], pts[min_i, 1])
                            annot.set_text(sel_labels[min_i])
                            annot.set_visible(True)
                            fig.canvas.draw_idle()
                    elif _hover_idx[0] >= 0:
                        _hover_idx[0] = -1
                        annot.set_visible(False)
                        fig.canvas.draw_idle()

                _current_indices = [[]]
                _shift_held = [False]

                def _update_selection(indices):
                    _current_indices[0] = indices
                    if indices:
                        sel_plot.set_data(pts[indices, 0], pts[indices, 1])
                    else:
                        sel_plot.set_data([], [])
                    fig.canvas.draw_idle()
                    layer_sel = defaultdict(list)
                    for i in indices:
                        lyr, fid = sel_fids[i]
                        layer_sel[id(lyr)].append(fid)
                    all_layers = {id(lyr): lyr for lyr, _ in sel_fids}
                    for lid, lyr in all_layers.items():
                        lyr.selectByIds(layer_sel[lid]) if lid in layer_sel else lyr.removeSelection()

                def _on_lasso(verts):
                    pts_axes = ax.transAxes.inverted().transform(
                        ax.transData.transform(pts))
                    new_idx = np.where(Path(verts).contains_points(pts_axes))[0].tolist()
                    if _shift_held[0]:
                        combined = list(set(_current_indices[0]) | set(new_idx))
                    else:
                        combined = new_idx
                    _update_selection(combined)

                _press_xy = [None]

                def _on_press(event):
                    if event.inaxes == ax and event.button == 1:
                        _press_xy[0] = (event.x, event.y)

                def _on_release(event):
                    if event.button != 1 or _press_xy[0] is None or event.xdata is None:
                        _press_xy[0] = None
                        return
                    moved = (abs(event.x - _press_xy[0][0]) > 5 or
                             abs(event.y - _press_xy[0][1]) > 5)
                    _press_xy[0] = None
                    if moved or event.inaxes != ax:
                        return
                    dists = np.hypot(pts[:, 0] - event.xdata, pts[:, 1] - event.ydata)
                    min_i = int(np.argmin(dists))
                    if dists[min_i] < 0.1:
                        if event.key == 'shift':
                            combined = list(set(_current_indices[0]) | {min_i})
                        else:
                            combined = [min_i]
                        _update_selection(combined)

                def _on_key_press(event):
                    if event.key == 'escape':
                        _update_selection([])
                    elif event.key == 'shift':
                        _shift_held[0] = True

                def _on_key_release(event):
                    if event.key == 'shift':
                        _shift_held[0] = False

                self._stereonet_lasso = _BoundedLassoSelector(ax, _on_lasso)
                fig.canvas.mpl_connect('button_press_event', _on_press)
                fig.canvas.mpl_connect('button_release_event', _on_release)
                fig.canvas.mpl_connect('key_press_event', _on_key_press)
                fig.canvas.mpl_connect('key_release_event', _on_key_release)
                fig.canvas.mpl_connect('motion_notify_event', _on_hover)

            ax.set_title(layer.name() + " [# " + str(len(iter)) + "]", pad=24)
            plt.show()

        else:
            self.iface.messageBar().pushMessage("No data selected, or no structural data found: first select a layer with structural info, then select the points that you wish to plot", level=Qgis.Warning, duration=5)
        
    def fieldExists(self, layer, fieldnames):
        """Backward-compatible wrapper for older internal calls."""
        return self._field_exists(layer, fieldnames)