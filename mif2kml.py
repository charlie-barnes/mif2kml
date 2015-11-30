#!/usr/bin/env python
#-*- coding: utf-8 -*-

#mif2kml - 
#This application is free software; you can redistribute
#it and/or modify it under the terms of the GNU General Public License
#defined in the COPYING file

#2008-2010 Charlie Barnes.

import sys
import os
import gtk
import gobject
import time
import sys
from math import pi, cos, tan, sin, pow, sqrt, floor, degrees, radians, atan2
import zipfile
from xml.sax.saxutils import escape
import csv
from optparse import OptionParser

class mif2kmlActions():
    def __init__(self):
    
        self.cancel = False
        
        if len(sys.argv) > 1:
        
            parser = OptionParser()
            parser.add_option("--mif", dest="mif",
                              help="MIF file", metavar="FILE")
            parser.add_option("--mid", dest="mid",
                              help="MID file", metavar="FILE")
            parser.add_option("--label", dest="label", metavar="INTEGER",
                              help="column number of the MID data to use as a label for the entities")
            parser.add_option("--sort", dest="sort", metavar="INTEGER",
                              help="column number of the MID data to sort by")
            parser.add_option("--directory", dest="output", metavar="DIR",
                              help="where the KML will be saved")
            parser.add_option("--kmz", dest="kmz", metavar="BOOLEAN",
                              help="compress the KML")

            (options, args) = parser.parse_args()
            
            self.process_file(mif=options.mif,
                              mid=options.mid,
                              kmz=options.kmz,
                              label=bool(options.label),
                              progress=False,
                              output=options.output,
                             )
            exit()

        #Load the widget tree
        builder = ""
        self.builder = gtk.Builder()
        self.builder.add_from_string(builder, len(builder))
        self.builder.add_from_file("ui.xml")

        signals = {
                   "mainQuit":self.main_quit,
                   "convert":self.convert,
                   "processFile":self.process_file,
                   "showAboutDialog":self.show_about_dialog,
                   "cancelConvert":self.cancel_convert,
                   "changeLabel":self.change_label,
                  }
        self.builder.connect_signals(signals)

        #Setup the main window
        self.main_window = self.builder.get_object("window1")
                
        self.main_window.show()
        self.create_filechooser_buttons()
        
    def create_filechooser_buttons(self):
        try:
            self.builder.get_object("eventbox1").get_child().destroy()
        except AttributeError:
            pass
        finally:
            filechooserbutton = gtk.FileChooserButton('Select a File')
            filechooserbutton.connect("file-set", self.allow_convert)
            self.builder.get_object("eventbox1").add(filechooserbutton)     
            filechooserbutton.show()
            
        try:
            self.builder.get_object("eventbox2").get_child().destroy()
        except AttributeError:
            pass
        finally:
            filechooserbutton = gtk.FileChooserButton('Select a File')
            filechooserbutton.connect("file-set", self.allow_spin)
            filechooserbutton.set_sensitive(False)
            self.builder.get_object("eventbox2").add(filechooserbutton)     
            filechooserbutton.show()
            
        try:
            self.builder.get_object("eventbox3").get_child().destroy()
        except AttributeError:
            pass
        finally:
            filechooserbutton = gtk.FileChooserButton('Select a Folder')
            filechooserbutton.set_action(gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER)
            filechooserbutton.set_sensitive(False)
            self.builder.get_object("eventbox3").add(filechooserbutton)
            filechooserbutton.show()
            
        self.builder.get_object("button1").set_sensitive(False)
        self.builder.get_object("hbox1").set_sensitive(False)
        self.builder.get_object("checkbutton1").set_sensitive(False)
        self.builder.get_object("checkbutton2").set_sensitive(False)
        self.builder.get_object("checkbutton1").set_active(False)
        self.builder.get_object("checkbutton2").set_active(False)
        self.builder.get_object("spinbutton1").set_value(0)
        self.builder.get_object("label7").set_text("")
        self.builder.get_object("eventbox3").get_child().set_filename(os.getenv('USERPROFILE') or os.getenv('HOME'))
            
    def change_label(self, widget):

        widget = self.builder.get_object("eventbox2").get_child()
        
        mid = widget.get_filename()
        
        if mid is not None:
            csvfile = open(mid)
            dialect = csv.Sniffer().sniff(csvfile.read(1024))
            csvfile.seek(0)
            mid_reader = csv.reader(csvfile, dialect)

            for row in mid_reader:
                self.builder.get_object("label7").set_text(row[int(self.builder.get_object("spinbutton1").get_value())])
                break
        
    def allow_spin(self, widget):
        self.builder.get_object("hbox1").set_sensitive(True)
        self.builder.get_object("checkbutton2").set_sensitive(True)

        mid = widget.get_filename()
        
        if mid is not None:
            csvfile = open(mid)
            dialect = csv.Sniffer().sniff(csvfile.read(1024))
            csvfile.seek(0)
            mid_reader = csv.reader(csvfile, dialect)

            for row in mid_reader:
                self.builder.get_object("spinbutton1").set_range(0, len(row)-1)
                break
                
        self.change_label(widget)
            
    def allow_convert(self, widget):
        self.builder.get_object("button1").set_sensitive(True)
        self.builder.get_object("eventbox2").get_child().set_sensitive(True)
        self.builder.get_object("eventbox3").get_child().set_sensitive(True)
        self.builder.get_object("eventbox2").get_child().set_current_folder(os.path.dirname(widget.get_filename()))
        self.builder.get_object("eventbox3").get_child().set_current_folder(os.path.dirname(widget.get_filename()))
        self.builder.get_object("checkbutton1").set_sensitive(True)

    def convert(self, widget):
        self.builder.get_object("progressbar1").show()
        self.builder.get_object("table1").set_sensitive(False)
        self.builder.get_object("button1").hide()
        self.builder.get_object("button6").show()
        self.builder.get_object("button2").set_sensitive(False)

        self.process_file(mif=self.builder.get_object("eventbox1").get_child().get_filename(), mid=self.builder.get_object("eventbox2").get_child().get_filename(), kmz=self.builder.get_object("checkbutton1").get_active(), label=self.builder.get_object("spinbutton1").get_value_as_int(),  sort=self.builder.get_object("checkbutton2").get_active(), progress=self.builder.get_object("progressbar1"), output=self.builder.get_object("eventbox3").get_child().get_filename())

        self.cancel = False
        self.builder.get_object("progressbar1").hide()
        self.builder.get_object("table1").set_sensitive(True)
        self.builder.get_object("button1").show()
        self.builder.get_object("button6").hide()
        self.builder.get_object("button2").set_sensitive(True)
        self.create_filechooser_buttons()
           
    def convertOSGB36toWGS84(self, lon, lat):

        lat = radians(lat)
        lon = radians(lon)

        a = 6377563.396
        b = 6356256.910

        sinPhi = sin(lat)
        cosPhi = cos(lat)
        sinLambda = sin(lon)
        cosLambda = cos(lon)
        H = 0

        eSq = (a*a - b*b) / (a*a)
        nu = a / sqrt(1 - eSq*sinPhi*sinPhi)

        x1 = (nu+H) * cosPhi * cosLambda
        y1 = (nu+H) * cosPhi * sinLambda
        z1 = ((1-eSq)*nu + H) * sinPhi

        # -- apply helmert transform using appropriate params
        tx = 446.448
        ty = -125.157
        tz = 542.060
        rx = 0.1502/3600 * pi/180  # normalise seconds to radians
        ry = 0.2470/3600 * pi/180
        rz = 0.8421/3600 * pi/180
        s1 = -20.4894/1e6 + 1              # normalise ppm to (s+1)

        # apply transform
        x2 = tx + x1*s1 - y1*rz + z1*ry
        y2 = ty + x1*rz + y1*s1 - z1*rx
        z2 = tz - x1*ry + y1*rx + z1*s1

        # -- convert cartesian to polar coordinates (using ellipse 2)
        a = 6378137
        b = 6356752.3142
        precision = 4 / a  # results accurate to around 4 metres

        eSq = (a*a - b*b) / (a*a)
        p = sqrt(x2*x2 + y2*y2)
        phi = atan2(z2, p*(1-eSq))
        phiP = 2*pi
        while (abs(phi-phiP) > precision):
            nu = a / sqrt(1 - eSq*sin(phi)*sin(phi))
            phiP = phi
            phi = atan2(z2 + eSq*nu*sin(phi), p)

        lam = atan2(y2, x2)
        H = p/cos(phi) - nu

        return degrees(lam), degrees(phi)
    
    def osgb_to_lonlat (self, east, north):
		# Airy 1830 major & minor semi-axes
	    a = 6377563.396 
	    b = 6356256.910 
	    # NatGrid scale factor on central meridian
	    F0 = 0.9996012717
	    # eccentricity squared                
	    e2 = 1 - (b**2)/(a**2)                          
	    n = (a-b)/(a+b)
	    n2 = n**2
	    n3 = n**3
	
	    lat=radians (49)
	    M=0
	    while (True):
		    lat = (north - -100000 - M)/(a*F0) + lat
		    Ma = (1 + n + (5/4)*n2 + (5/4)*n3) * (lat-radians (49))
		    Mb = (3*n + 3*n2 + (21/8)*n3) * sin (lat-radians (49)) * cos (lat+radians (49))
		    Mc = ((15/8)*n2 + (15/8)*n3) * sin (2*(lat-radians (49))) * cos (2*(lat+radians (49)))
		    Md = (35/24)*n3 * sin(3*(lat-radians (49))) * cos (3*(lat+radians (49)))
		    # meridional arc
		    M = b * F0 * (Ma - Mb + Mc - Md)
		    if (north - -100000 - M <= 0.00001):
			    # repeat until < 0.01mm
			    break
			
	    sinlat = sin(lat)
	    # transverse radius of curvature
	    nu = a*F0 / sqrt (1-e2*sinlat*sinlat)
	    # meridional radius of curvature     
	    rho = a * F0 * (1 - e2) / pow (1 - e2 * sinlat**2, 1.5)
	    eta2 = nu / rho - 1
	    tanlat = tan (lat)
	    tanlat2 = tanlat**2
	    tanlat4 = tanlat2**2
	    tanlat6 = tanlat4 * tanlat2
	    seclat = 1 / cos (lat)
	    nu3 = nu**3
	    nu5 = nu3 * nu**2
	    nu7 = nu5 * nu**2
	    VII = tanlat / (2*rho*nu)
	    VIII = tanlat / (24*rho*nu3) * (5+3*tanlat2+eta2-9*tanlat2*eta2)
	    IX = tanlat / (720*rho*nu5) * (61+90*tanlat2+45*tanlat4)
	    X = seclat / nu
	    XI = seclat / (6*nu3) * (nu/rho+2*tanlat2)
	    XII = seclat / (120*nu5) * (5+28*tanlat2+24*tanlat4)
	    XIIA = seclat / (5040*nu7) * (61+662*tanlat2+1320*tanlat4+720*tanlat6)
	    dE = east - 400000
	    lat = lat - VII*dE**2 + VIII*dE**4 - IX*dE**6
	    lon = radians (-2) + X*dE - XI*dE**3 + XII*dE**5 - XIIA*dE**7

	    return degrees(lon), degrees(lat)


    def process_file(self, mif, mid=None, kmz=False, label=0, sort=False, progress=False, output=None):

        if kmz == 'True' or kmz == True:
            kmz = True

        if mif is None:
            print "No MIF file to convert"
            exit()

        # hacky remove extension!
        if os.path.basename(mif)[-4:-3] == ".":    
            name = os.path.basename(mif)[:-4]
        else:
            name = os.path.basename(mif)

        mifsource = mif
        
        try:    
            mif = open(mifsource, 'r')
        except IOError:
            print "Can't find that file"
            exit()

        if output is None:
            directory = os.path.dirname(mifsource)
        else:
            directory = output
            
        entity_names = {}
        sort_names = {}

        if mid is not None:
            try:    
                csvfile = open(mid)
                dialect = csv.Sniffer().sniff(csvfile.read(1024))
                csvfile.seek(0)
                mid_reader = csv.reader(csvfile, dialect)

                linecount = 0
            
                for row in mid_reader:
                    linecount = linecount+1
                    try:
                        entity_names[linecount] = escape(row[label])
                    except:
                        entity_names[linecount] = escape(row[0])
            except IOError:
                print "MID file not found"
                mid = None

        try:
            kml = open(''.join([os.path.join(directory, name), ".kml"]), 'w')
        except:
            print "Unable to write to output directory"
            return -2
            
        totallinecount = 2
        for line in mif:
            totallinecount = totallinecount + 1.0
            
        mif.seek(0)

        in_data = False
        in_geom = False
        count = 0
        subcoun = 0
        linecount = 0
        coord = None


        entities = {}
         
        kml.write(''.join(["<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n<Folder>\n<name>", escape(name), "</name><open>1</open>"]))

        if progress == False:
            print "Converting..."

        for line in mif:
            if self.cancel:
                kml.close()
                os.remove(''.join([os.path.join(directory, name), ".kml"]))
                return -4
                
            linecount = linecount + 1.0

            if progress:
                if linecount%100==0:
                    progress.set_fraction(linecount/totallinecount)

                    while gtk.events_pending():
                       gtk.main_iteration()

            if line[0:8].upper() == "COORDSYS": # check the coord system
                if line[0:31].upper() == "COORDSYS EARTH PROJECTION 8, 79":
                    coord = "bng"
                elif line[0:32].upper() == "COORDSYS EARTH PROJECTION 1, 104":
                    coord = "wgs84"
                else:
                    return -1
            elif line[0:4].upper() == "DATA": # we're into the data
                in_data = True
            elif in_data:      
                if (line[0:2] == "  ") and (line[2:3].isdigit()): # we've come to a "sub polygon"
                    subcount = subcount + 1
                    entities[entity_name].append([])
                elif (line[0:6].upper() == "REGION"): # we've found a new region entity
                
                    subcount = 0
                    in_geom = True
                    
                    count = count+1
                    
                    if label >= 0:
                        entity_name = entity_names[count]
                    else:
                        entity_name = ''.join(["Region ", str(count)])
                    
                    entities[entity_name] = [[]]
                    
                    mif.next() #ignore the next line (the polygon vertices count)
                    
                elif in_geom and (line[0:1] != " "):
                    if coord == "bng":
                        x = int(float(line.split(' ')[0].strip()))
                        y = int(float(line.split(' ')[1].strip()))
               
                        lon, lat = self.osgb_to_lonlat(x, y) #bng to osgb36
                        lon, lat = self.convertOSGB36toWGS84(lon, lat) #osgb36 to wgs84
                    elif coord == "wgs84":
                        lon = line.split(' ')[0].strip()
                        lat = line.split(' ')[1].strip()
                        
                    entities[entity_name][subcount].append([lon, lat])


        keys = entities.keys()
        if sort:
            keys.sort()

        for entity in keys:
            kml.write(''.join(["\n<Placemark>\n<name>", entity , "</name>\n<MultiGeometry>\n<Polygon>\n<outerBoundaryIs>\n<LinearRing>\n<coordinates>\n"]))
            
            for subentities in entities[entity]:
                for points in subentities:
                    kml.write(''.join([str(points[0]), ",", str(points[1]), ",",  "0\n"]))

                kml.write("\n</coordinates>\n</LinearRing>\n</outerBoundaryIs>\n</Polygon>\n<Polygon>\n<outerBoundaryIs>\n<LinearRing>\n<coordinates>\n")

            kml.write("</coordinates>\n</LinearRing>\n</outerBoundaryIs>\n</Polygon>\n</MultiGeometry>\n</Placemark>\n")

        kml.write("</Folder>\n</kml>")
        kml.close()

        if kmz == True:
            if progress == False:
                print "Compressing..."
            kmz = zipfile.ZipFile(''.join([os.path.join(directory, name), ".kmz"]), "w")
            kmz.write(''.join([os.path.join(directory, name), ".kml"]), ''.join([name, ".kml"]), zipfile.ZIP_DEFLATED)
            os.remove(''.join([os.path.join(directory, name), ".kml"]))

    def main_quit(self, widget, var=None):
        gtk.main_quit()

    def cancel_convert(self, widget):
       self.cancel = True

    def show_about_dialog(self, widget):
       about=gtk.AboutDialog()
       about.set_name("mif2kml")
       about.set_copyright("2010 Charlie Barnes")
       about.set_authors(["Charlie Barnes <charlie@cucaera.co.uk>"])
       about.set_license("mif2kml is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the Licence, or (at your option) any later version.\n\nmif2kml is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.\n\nYou should have received a copy of the GNU General Public License along with mif2kml; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA")
       about.set_wrap_license(True)
       about.set_website("http://cucaera.co.uk/software/mif2kml/")
       about.set_transient_for(self.builder.get_object("window1"))
       result=about.run()
       about.destroy()

if __name__ == '__main__':
    mif2kmlActions()
    gtk.main()
    
