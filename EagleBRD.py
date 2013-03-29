# * File: EagleBRD.py **********************************************************
# *                                                                            *
# * copyright (c) 2011 Franz Riedmueller <riedmueller.f@gmail.com>             *
# *                                                                            *
# * This program is free software: you can redistribute it and/or modify       *
# * it under the terms of the GNU General Public License as published by       *
# * the Free Software Foundation, either version 3 of the License, or          *
# * (at your option) any later version.                                        *
# *                                                                            *
# * This program is distributed in the hope that it will be useful,            *
# * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
# * GNU General Public License for more details.                               *
# *                                                                            *
# * You should have received a copy of the GNU General Public License          *
# * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
# *                                                                            *
# * Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen   *
# * der GNU General Public License, wie von der Free Software Foundation,      *
# * Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren              *
# * veroeffentlichten Version, weiterverbreiten und/oder modifizieren.         *
# *                                                                            *
# * Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber    *
# * OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite        *
# * Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FUER EINEN BESTIMMTEN ZWECK
# * Siehe die GNU General Public License fuer weitere Details.                 *
# *                                                                            *
# * Sie sollten eine Kopie der GNU General Public License zusammen mit diesem  *
# * Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>. *
# *                                                                            *
# ******************************************************************************

import FreeCAD
import Part
import Mesh
import math
import xml.dom.minidom
import string
import re
from FreeCAD import Base

import FreeCADGui, __builtin__
import os
import ImportGui

from ConfigParser import SafeConfigParser

# --- global define values -----------------------------------------------------
NOBOARD=0
PLANEGEN=1
CADFILE=2
EXTRACTED=3
PH_NONE=0
PH_BALL=1
PH_CYL=2

# --- std. values if no .cfg was used - Here you can set global default settings
default_docname="NEWDOC"
default_boardgentype=PLANEGEN
default_thickness=1.6
default_show_notfound=PH_BALL
default_bgen_border=10.
default_libpath=FreeCAD.getHomePath()+"Mod/EagleBRD/lib/"
default_pack_excludelist=["bla"]
default_signalgen=True
default_diagnostic=True
default_diag_path = FreeCAD.getHomePath()+"Mod/EagleBRD/"
default_diag_filename = "diagnostic.txt"
default_signal_pcb_copper = 0.02
default_signal_wirelimit = 4000

# --- global variables ---------------------------------------------------------
docname = default_docname
boardgentype = default_boardgentype
pcb_thickness = default_thickness
show_notfound = default_show_notfound
bgen_border = default_bgen_border
libpath = default_libpath
pack_excludelist = []
signalgen = default_signalgen
boardCADfilePath = ""
#boardCADfilePath = libpath + "board.stl"
libFilePath = libpath+"translib.xml"
USR_libFilePath = libpath+"USR_translib.xml"
DIAGNOSTIC = default_diagnostic
diag_path = default_diag_path
diag_filename = default_diag_filename
lib_pack = []
lib_file = []
lib_desc = []
lib_movx = []
lib_movy = []
lib_movz = []
lib_rotx = []
lib_roty = []
lib_rotz = []
signal_pcb_copper = default_signal_pcb_copper
signal_wirelimit = default_signal_wirelimit

# currently not used
POINT_TOLERANCE = 0.005
LAYER_3D_NUMBER = "57"
#inputFilePath = 'D:/bin/FreeCAD/Mod/MyScripts/test.brd'
#dom = xml.dom.minidom.parse(inputFilePath)
#domlib = xml.dom.minidom.parse(libFilePath)


# ******************************************************************************
# This function was called after opening a Eagle .brd file with FreeCAD
# == main()
def open(filename):
    """called when freecad opens an BRD file"""
    #docname1 = os.path.splitext(os.path.basename(filename))[0]
    #doc = FreeCAD.newDocument(docname1)
    #FreeCAD.setActiveDocument(docname1)
    #message='Started with opening of "'+filename+'" file\n'
    FreeCAD.Console.PrintMessage(filename)
    configsystem(filename)
    dom = xml.dom.minidom.parse(filename)
    placeparts(dom)


# ******************************************************************************
# This function was called after importing a Eagle .brd file with FreeCAD
# == main()
def insert(filename,docname1):
    """called when freecad imports an BRD file"""
    #FreeCAD.setActiveDocument(docname1)
    #doc=FreeCAD.getDocument(docname1)
    #FreeCAD.Console.PrintMessage('Started import of "'+filename+'" file')
    FreeCAD.Console.PrintMessage(filename)
    configsystem(filename)
    dom = xml.dom.minidom.parse(filename)
    placeparts(dom)


# ******************************************************************************
# This function does the standard configuration and reads out an optinal .cfg
# file for overriding the standard values
def configsystem(fname):
  global pcb_thickness
  global docname
  global boardgentype
  global bgen_border
  global pack_excludelist
   
  docname = string.replace(os.path.splitext(os.path.basename(fname))[0],'.','_')
  docname = string.replace(docname,'-','_')
  #print os.path.dirname(fname)
  if DIAGNOSTIC==True:
    diagfile=__builtin__.open(diag_path+diag_filename, "w")
    diagfile.writelines("Working path: "+os.path.dirname(fname)+"\n")
    diagfile.writelines("Configfile: "+os.path.dirname(fname)+"/"+os.path.splitext(os.path.basename(fname))[0]+".cfg\n")
    diagfile.close()
  
  cfgfile=os.path.dirname(fname)+"/"+os.path.splitext(os.path.basename(fname))[0]+".cfg"
  if os.path.isfile(cfgfile):
    if DIAGNOSTIC==True:
      diagfile=__builtin__.open(diag_path+diag_filename, "a")
      diagfile.writelines("Found: Config File\n")
      diagfile.close()

    parser = SafeConfigParser()
    parser.read(cfgfile)
       
    if parser.has_section('BOARD'):
      if DIAGNOSTIC==True:
        diagfile=__builtin__.open(diag_path+diag_filename, "a")
        diagfile.writelines("Found: Section BOARD in .cfg file\n")
        diagfile.close()
      if parser.has_option('BOARD','BOARD_THICKNESS'):
        pcb_thickness = parser.getfloat('BOARD','BOARD_THICKNESS')
      else:
        pcb_thickness = default_thickness
        
      if parser.has_option('BOARD','BOARD_GENTYPE'):
        boardgentype = parser.getint('BOARD','BOARD_GENTYPE')
        if (boardgentype>EXTRACTED):
          boardgentype=NOBOARD
    if parser.has_section('LIBRARY'):
      if parser.has_option('LIBRARY','LIB_PATH'):
        libpath = parser.getstring('LIBRARY','LIB_PATH')
      else:
        libpath = default_libpath
    if parser.has_section('PARTS'):
      if parser.has_option('PARTS','EXCLUDELIST'):
        list = (parser.get('PARTS','EXCLUDELIST')).split(",")
        for i in range(0,len(list)):
          list[i] = re.sub(r'\s', '', list[i])
        for ele in list:
          if ele=="":
            list.remove(ele)
          else:
            pack_excludelist.append(ele)
            
        if DIAGNOSTIC==True:
          diagfile=__builtin__.open(diag_path+diag_filename, "a")
          diagfile.writelines("Found: Option EXCLUDELIST in Section PARTS in .cfg\n")
          for xy in pack_excludelist:
            diagfile.writelines("pack_excludelist Element:"+xy+"\n")
            #diagfile.writelines("EXCLUDELIST Element:"+xy+"\n")
          diagfile.close()
      else:
        pack_excludelist=default_pack_excludelist
    if parser.has_section('BOARDGEN'):
      if parser.has_option('BOARDGEN','BGEN_BORDER'):
        bgen_border = parser.getfloat('BOARDGEN','BGEN_BORDER')
      else:
        bgen_border = default_bgen_border
    if parser.has_section('SIGNALS'):
      if parser.has_option('SIGNALS','SIG_ON'):
        signalgen = True
      else:
        signalgen = False

  else:
    if DIAGNOSTIC==True:
      diagfile=__builtin__.open(diag_path+diag_filename, "a")
      diagfile.writelines("Config file not found. Using standard Values.\n")
      diagfile.close()


# ******************************************************************************
# This function converts euler angle format to FreeCAD Quaternion format
def eulerToQuaternion(yaw=0.0, pitch=0.0, roll=0.0):
  yaw = math.radians(yaw)
  pitch = math.radians(pitch)
  roll = math.radians(roll)
  
  c1 = math.cos(yaw / 2.)
  c2 = math.cos(pitch / 2.)
  c3 = math.cos(roll / 2.)
  s1 = math.sin(yaw / 2.)
  s2 = math.sin(pitch / 2.)
  s3 = math.sin(roll / 2.)
  
  q1 = c1 * c2 * s3 - s1 * s2 * c3
  q2 = c1 * s2 * c3 + s1 * c2 * s3
  q3 = s1 * c2 * c3 - c1 * s2 * s3
  q4 = c1 * c2 * c3 + s1 * s2 * s3

  return FreeCAD.Rotation(q1, q2, q3, q4)


# ******************************************************************************
# This function parses the library-translation-file an sets the depending 
# variables an arrays
def parselibfile():
  if os.path.isfile(USR_libFilePath):
    USR_domlib = xml.dom.minidom.parse(USR_libFilePath)
    USR_lib1 = USR_domlib.getElementsByTagName("translib")[0]
    translations = USR_lib1.getElementsByTagName("trans")
    for t in translations:
      lib_pack.append(t.getAttribute("package"))
      lib_file.append(t.getAttribute("file"))
      lib_desc.append(t.getAttribute("desc"))
      lib_movx.append(t.getAttribute("movx"))
      lib_movy.append(t.getAttribute("movy"))
      lib_movz.append(t.getAttribute("movz"))
      lib_rotx.append(t.getAttribute("rotx"))
      lib_roty.append(t.getAttribute("roty"))
      lib_rotz.append(t.getAttribute("rotz"))
  #-->standard translations
  domlib = xml.dom.minidom.parse(libFilePath)
  lib1 = domlib.getElementsByTagName("translib")[0]
  translations = lib1.getElementsByTagName("trans")
  for t in translations:
      lib_pack.append(t.getAttribute("package"))
      lib_file.append(t.getAttribute("file"))
      lib_desc.append(t.getAttribute("desc"))
      lib_movx.append(t.getAttribute("movx"))
      lib_movy.append(t.getAttribute("movy"))
      lib_movz.append(t.getAttribute("movz"))
      lib_rotx.append(t.getAttribute("rotx"))
      lib_roty.append(t.getAttribute("roty"))
      lib_rotz.append(t.getAttribute("rotz"))


# ******************************************************************************
# This function converts ARC values  from Eagle to Values for FreeCAD
# input:
# ax1/ay1 == startingpoint, ax2/ay2 == endpoint, aangle=included angle in degree
# output:
# centerx/y == coordinates of circle, radius=radius of circle, 
# startangle/endangle are the angles of the arc in degree
#
def convertarcvalues(ax1,ay1,ax2,ay2,aangle):
  dist1 = math.sqrt(math.pow((ax1-ax2),2)+math.pow((ay1-ay2),2))
  # Radius of ARC
  radius = dist1/math.sqrt(2.0*(1.0-math.cos(math.radians(aangle))))
  # Middlepoint between the two given points
  xm=(ax1+ax2)/2.
  ym=(ay1+ay2)/2.
  # slope through middlepoint and center of arc
  m=(ax1-ax2)/(ay2-ay1)
  # distance from center of arc to middlepoint
  a=math.sqrt(math.pow(radius,2)-math.pow((dist1/2.0),2))
  # center
  dx=(a/math.sqrt(math.pow(m,2)+1.0))
  dy=((m*a)/math.sqrt(math.pow(m,2)+1.0))
  if ((aangle<0) and (((ay1<ay2) and (ax1<ax2)) or ((ay2>ay1) and (ax2<ax1)))) or ((aangle>0) and (((ay1>ay2) and (ax1>ax2)) or ((ay2<ay1) and (ax2>ax1)))):
    centery = ym+dy
    centerx = xm+dx
  else:
    centery = ym-dy
    centerx = xm-dx
  startangle = math.atan2((ay1-centery),(ax1-centerx))*180.0/math.pi
  endangle = startangle + aangle
  if (startangle < 0.0) or (endangle<0.0):
    startangle = startangle + 360.0
    endangle = endangle + 360.0
  if startangle > endangle:
    mkr = startangle
    startangle = endangle
    endangle = mkr
  return centerx, centery, radius, startangle, endangle
  #return centerx, centery, radius, startangle, endangle, xm, ym, dx, dy, m, a


# ******************************************************************************
# Func converts ARC values from Eagle to Values for FreeCAD ARC (threepoints)
# input:
# ax1/ay1 == startingpoint, ax2/ay2 == endpoint, aangle=included angle in degree
# output:
# ax1/ax2/ay1/ay2 == coordinates of the start and endpoint of ARC,
# a_middle_x/_y == Middlepoint on the ARC
#
def conv_arc_to_threepoint(ax1,ay1,ax2,ay2,aangle):
  Cx, Cy, rad, anglestart, angleend = convertarcvalues(ax1,ay1,ax2,ay2,aangle)
  a_middle_x = float(rad) * math.cos( math.radians((float(angleend-anglestart))/2. + float(anglestart)) )
  a_middle_y = float(rad) * math.sin( math.radians((float(angleend-anglestart))/2. + float(anglestart)) )
  a_middle_x = a_middle_x + Cx
  a_middle_y = a_middle_y + Cy
  return ax1, ay1, a_middle_x, a_middle_y, ax2, ay2


# ******************************************************************************
# Generates signals and VIA's out of *.brd copper layers and Group it into 
# signalnames and layers
def showsignals(bfile1):
  global signal_pcb_copper
  global signal_wirelimit
  lay1_shapes = []
  lay2_shapes = []
  lay15_shapes = []
  lay16_shapes =[]
  wirecnt = 0
  # Get the Eagle Board File signals
  path1 = bfile1.getElementsByTagName("eagle")[0]
  path2 = path1.getElementsByTagName("drawing")[0]
  path3 = path2.getElementsByTagName("board")[0]
  path4 = path3.getElementsByTagName("signals")[0]
  signals = path4.getElementsByTagName("signal")
  
  FreeCAD.activeDocument().addObject("App::DocumentObjectGroup","SIGNALS")
  
  for sg in signals:
    while len(lay1_shapes)>0:
      lay1_shapes.pop()
    while len(lay2_shapes)>0:
      lay2_shapes.pop()
    while len(lay15_shapes)>0:
      lay15_shapes.pop()
    while len(lay16_shapes)>0:
      lay16_shapes.pop()
    sgname = sg.getAttribute("name")
    sgwires = sg.getElementsByTagName("wire")
    layer1_not_empty = 0
    layer2_not_empty = 0
    layer15_not_empty = 0
    layer16_not_empty = 0
    for sgw in sgwires:
      if wirecnt < signal_wirelimit:
        if ((sgw.getAttribute("layer")=="1") or (sgw.getAttribute("layer")=="2") or (sgw.getAttribute("layer")=="15") or (sgw.getAttribute("layer")=="16")):
          wirecnt = wirecnt+1
          sgx1 = float(sgw.getAttribute("x1"))
          sgy1 = float(sgw.getAttribute("y1"))
          sgx2 = float(sgw.getAttribute("x2"))
          sgy2 = float(sgw.getAttribute("y2"))
          sgwidth = float(sgw.getAttribute("width"))
          if (sgw.hasAttribute("curve")!=True):
            sglen = math.sqrt(math.pow((sgx1-sgx2),2)+math.pow((sgy1-sgy2),2))
            box = Part.makeBox(sglen,sgwidth,signal_pcb_copper,FreeCAD.Vector(0,-(sgwidth/2),0))
            cyl1= Part.makeCylinder((sgwidth/2),signal_pcb_copper)
            cyl2= Part.makeCylinder((sgwidth/2),signal_pcb_copper,FreeCAD.Vector(sglen,0,0))
            wireseg = box
            wireseg = wireseg.fuse(cyl1)
            wireseg = wireseg.fuse(cyl2)
            segangle = math.atan2((sgy2-sgy1),(sgx2-sgx1))*180/math.pi
            if sgw.getAttribute("layer")=="1":
              layer1_not_empty = 1
              wireseg.Placement = FreeCAD.Placement(FreeCAD.Vector(sgx1, sgy1, (0.5*pcb_thickness)), FreeCAD.Rotation(eulerToQuaternion(segangle, 0, 0)))
              lay1_shapes.append(wireseg)
            elif sgw.getAttribute("layer")=="2":
              layer2_not_empty = 1
              wireseg.Placement = FreeCAD.Placement(FreeCAD.Vector(sgx1, sgy1, (pcb_thickness/6)), FreeCAD.Rotation(eulerToQuaternion(segangle, 0, 0)))
              lay2_shapes.append(wireseg)
            elif sgw.getAttribute("layer")=="15":
              layer15_not_empty = 1
              wireseg.Placement = FreeCAD.Placement(FreeCAD.Vector(sgx1, sgy1, -(pcb_thickness/6)), FreeCAD.Rotation(eulerToQuaternion(segangle, 0, 0)))
              lay15_shapes.append(wireseg)
            elif sgw.getAttribute("layer")=="16":
              layer16_not_empty = 1
              wireseg.Placement = FreeCAD.Placement(FreeCAD.Vector(sgx1, sgy1, (-0.5*pcb_thickness)-signal_pcb_copper), FreeCAD.Rotation(eulerToQuaternion(segangle, 0, 0)))
              lay16_shapes.append(wireseg)
          #elif (False==True):
          elif (sgw.hasAttribute("curve")==True):
            arccurve = float(sgw.getAttribute("curve"))
            cyl1 = Part.makeCylinder((sgwidth/2),signal_pcb_copper)
            cyl2 = Part.makeCylinder((sgwidth/2),signal_pcb_copper)
            arcx,arcy,arcrad,arcstart,arcstop = convertarcvalues(sgx1,sgy1,sgx2,sgy2,arccurve)
            if DIAGNOSTIC==True:
              diagfile=__builtin__.open(diag_path+diag_filename, "a")
              diagfile.writelines("ARC: x="+str(arcx)+"  y="+str(arcy)+"  rad="+str(arcrad)+"  start="+str(arcstart)+"  stop="+str(arcstop)+"\n")
              diagfile.close()
            if arccurve<0:
              arccurve=-arccurve
            arc1 = Part.makeCylinder(arcrad+(sgwidth/2),signal_pcb_copper,FreeCAD.Vector(0,0,0),FreeCAD.Vector(0,0,1),arccurve)
            arc1 = arc1.cut(Part.makeCylinder(arcrad-(sgwidth/2),signal_pcb_copper,FreeCAD.Vector(0,0,0),FreeCAD.Vector(0,0,1),arccurve))
            if sgw.getAttribute("layer")=="1":
              layer1_not_empty = 1
              arc1.Placement = FreeCAD.Placement(FreeCAD.Vector(arcx, arcy, (0.5*pcb_thickness)), FreeCAD.Rotation(eulerToQuaternion(arcstart, 0, 0)))
              cyl1.Placement = FreeCAD.Placement(FreeCAD.Vector(sgx1, sgy1, (0.5*pcb_thickness)), FreeCAD.Rotation(eulerToQuaternion(0, 0, 0)))
              cyl2.Placement = FreeCAD.Placement(FreeCAD.Vector(sgx2, sgy2, (0.5*pcb_thickness)), FreeCAD.Rotation(eulerToQuaternion(0, 0, 0)))
              arc1=arc1.fuse(cyl1)
              arc1=arc1.fuse(cyl2)
              lay1_shapes.append(arc1)
            elif sgw.getAttribute("layer")=="2":
              layer2_not_empty = 1
              arc1.Placement = FreeCAD.Placement(FreeCAD.Vector(arcx, arcy, (pcb_thickness/6)), FreeCAD.Rotation(eulerToQuaternion(arcstart, 0, 0)))
              cyl1.Placement = FreeCAD.Placement(FreeCAD.Vector(sgx1, sgy1, (pcb_thickness/6)), FreeCAD.Rotation(eulerToQuaternion(0, 0, 0)))
              cyl2.Placement = FreeCAD.Placement(FreeCAD.Vector(sgx2, sgy2, (pcb_thickness/6)), FreeCAD.Rotation(eulerToQuaternion(0, 0, 0)))
              arc1=arc1.fuse(cyl1)
              arc1=arc1.fuse(cyl2)
              lay2_shapes.append(arc1)
            elif sgw.getAttribute("layer")=="15":
              layer15_not_empty = 1
              arc1.Placement = FreeCAD.Placement(FreeCAD.Vector(arcx, arcy, -(pcb_thickness/6)), FreeCAD.Rotation(eulerToQuaternion(arcstart, 0, 0)))
              cyl1.Placement = FreeCAD.Placement(FreeCAD.Vector(sgx1, sgy1, -(pcb_thickness/6)), FreeCAD.Rotation(eulerToQuaternion(0, 0, 0)))
              cyl2.Placement = FreeCAD.Placement(FreeCAD.Vector(sgx2, sgy2, -(pcb_thickness/6)), FreeCAD.Rotation(eulerToQuaternion(0, 0, 0)))
              arc1=arc1.fuse(cyl1)
              arc1=arc1.fuse(cyl2)
              lay15_shapes.append(arc1)
            elif sgw.getAttribute("layer")=="16":
              layer16_not_empty = 1
              arc1.Placement = FreeCAD.Placement(FreeCAD.Vector(arcx, arcy, (-0.5*pcb_thickness)-signal_pcb_copper), FreeCAD.Rotation(eulerToQuaternion(arcstart, 0, 0)))
              cyl1.Placement = FreeCAD.Placement(FreeCAD.Vector(sgx1, sgy1, (-0.5*pcb_thickness)-signal_pcb_copper), FreeCAD.Rotation(eulerToQuaternion(0, 0, 0)))
              cyl2.Placement = FreeCAD.Placement(FreeCAD.Vector(sgx2, sgy2, (-0.5*pcb_thickness)-signal_pcb_copper), FreeCAD.Rotation(eulerToQuaternion(0, 0, 0)))
              arc1=arc1.fuse(cyl1)
              arc1=arc1.fuse(cyl2)
              lay16_shapes.append(arc1)

    sgvias = sg.getElementsByTagName("via")
    for sgw in sgvias:
      if sgw.getAttribute("extent")=="1-16":
        layer1_not_empty = 1
        layer16_not_empty = 1
        vx = float(sgw.getAttribute("x"))
        vy = float(sgw.getAttribute("y"))
        vdrill = float(sgw.getAttribute("drill"))
        restring = 0.2
        # Via ring auf top seite
        vcyl= Part.makeCylinder((vdrill/2)+restring,signal_pcb_copper,FreeCAD.Vector(vx,vy,(0.5*pcb_thickness)))
        vcyl= vcyl.cut(Part.makeCylinder((vdrill/2),signal_pcb_copper,FreeCAD.Vector(vx,vy,(0.5*pcb_thickness))))
        lay1_shapes.append(vcyl)
        # Durchkontaktierung (Huelse) --> top zugehoerigkeit
        vcyl= Part.makeCylinder((vdrill/2)+signal_pcb_copper,pcb_thickness,FreeCAD.Vector(vx,vy,-(0.5*pcb_thickness)))
        vcyl= vcyl.cut(Part.makeCylinder((vdrill/2),pcb_thickness,FreeCAD.Vector(vx,vy,-(0.5*pcb_thickness))))
        lay1_shapes.append(vcyl)
        # Via ring auf bottom seite
        vcyl= Part.makeCylinder((vdrill/2)+restring,signal_pcb_copper,FreeCAD.Vector(vx,vy,(-0.5*pcb_thickness)-signal_pcb_copper))
        vcyl= vcyl.cut(Part.makeCylinder((vdrill/2),signal_pcb_copper,FreeCAD.Vector(vx,vy,(-0.5*pcb_thickness)-signal_pcb_copper)))
        lay16_shapes.append(vcyl)
    
    if(layer1_not_empty==1):
      comp=Part.Compound(lay1_shapes)
      Part.show(comp)
      FreeCAD.getDocument(FreeCAD.activeDocument().Name).ActiveObject.Label=sgname+"_TSIG"
      FreeCADGui.getDocument(docname).ActiveObject.ShapeColor=(0.,(85./255.),0.)
      FreeCAD.activeDocument().getObject("SIGNALS").addObject(FreeCAD.getDocument(FreeCAD.activeDocument().Name).ActiveObject)
    if(layer2_not_empty==1):
      comp=Part.Compound(lay2_shapes)
      Part.show(comp)
      FreeCAD.getDocument(FreeCAD.activeDocument().Name).ActiveObject.Label=sgname+"_IN1SIG"
      FreeCADGui.getDocument(docname).ActiveObject.ShapeColor=((218./255.),(138./255.),(103./255.))
      FreeCAD.activeDocument().getObject("SIGNALS").addObject(FreeCAD.getDocument(FreeCAD.activeDocument().Name).ActiveObject)
    if(layer15_not_empty==1):
      comp=Part.Compound(lay15_shapes)
      Part.show(comp)
      FreeCAD.getDocument(FreeCAD.activeDocument().Name).ActiveObject.Label=sgname+"_IN2SIG"
      FreeCADGui.getDocument(docname).ActiveObject.ShapeColor=((218./255.),(138./255.),(103./255.))
      FreeCAD.activeDocument().getObject("SIGNALS").addObject(FreeCAD.getDocument(FreeCAD.activeDocument().Name).ActiveObject)
    if(layer16_not_empty==1):
      comp=Part.Compound(lay16_shapes)
      Part.show(comp)
      FreeCAD.getDocument(FreeCAD.activeDocument().Name).ActiveObject.Label=sgname+"_BSIG"
      FreeCADGui.getDocument(docname).ActiveObject.ShapeColor=(0.,(85./255.),0.)
      FreeCAD.activeDocument().getObject("SIGNALS").addObject(FreeCAD.getDocument(FreeCAD.activeDocument().Name).ActiveObject)


# ******************************************************************************
# Opens a secondary FreeCAD FCStd file and copies all visible objects to the 
# current FreeCAD document (should be a kind of import functionality)
def insertfcstd(fcstdpath, dname):
  mydoc=FreeCAD.activeDocument()
  mydocname=FreeCAD.activeDocument().Name
  opendoc=FreeCAD.openDocument(fcstdpath)
  opendocname=FreeCAD.activeDocument().Name
  for i in opendoc.Objects:
    if (i.ViewObject.Visibility==True):
      cpy = mydoc.copyObject(i)
      cpy.ViewObject.DiffuseColor = i.ViewObject.DiffuseColor
      #cpy.ViewObject.ShapeColor = i.ViewObject.ShapeColor
      #cpy.ViewObject.ShapeMaterial = i.ViewObject.ShapeMaterial
  #FreeCAD.setActiveDocument(mydocname)
  #FreeCAD.getDocument(mydocname).ActiveObject.Label = newlabel
  FreeCAD.closeDocument(opendocname)
  FreeCAD.setActiveDocument(mydocname)


# ******************************************************************************
# Extracts a Part out of boardfile lib
# Currently not used / activated
def partextract(brdfile, plabel, pname, fusetol=0.0):
  extrusions = []
  edges = []

  while len(edges)>0:
    edges.pop()
  while len(extrusions)>0:
    extrusions.pop()
  
  path1 = brdfile.getElementsByTagName("eagle")[0]
  path2 = path1.getElementsByTagName("drawing")[0]
  path3 = path2.getElementsByTagName("board")[0]
  path4 = path3.getElementsByTagName("libraries")[0]
  libs = path4.getElementsByTagName("library")
  
  for L in libs:
    packages = L.getElementsByTagName("packages")[0]
    packs = packages.getElementsByTagName("package")
    
    for P in packs:
      if P.getAttribute("name")==pname:
        # TODO circles und nicht zusammenhaengende Strukturen
        # TODO geometrien mit gleicher hoehe aber zusammengesetzt (pins)
        wires = P.getElementsByTagName("wire")
        #edges=[]
        for W in wires:
          if W.getAttribute("layer")==LAYER_3D_NUMBER:
            # print W.getAttribute("width")
            # TODO: Alle Punkte auslesen und auf geschlossenes Polygon pruefen
            # TODO: ARC import
            # TODO: Circle import
            x1 = float(W.getAttribute("x1"))
            y1 = float(W.getAttribute("y1"))
            x2 = float(W.getAttribute("x2"))
            y2 = float(W.getAttribute("y2"))
            h  = 1000. * float(W.getAttribute("width"))
            vkt1=Base.Vector(x1,y1,0)
            vkt2=Base.Vector(x2,y2,0)
            L1 = Part.Line(vkt1,vkt2)
            #print L1
            edges.append(L1)

        shape = Part.Shape(edges)
        for v in shape.Vertexes: v.setTolerance(fusetol)

        W1 = Part.Wire(shape.Edges)
        F1 = Part.Face(W1)

        E1 = F1.extrude(Base.Vector(0,0,h))
        #Part.show(E1)
        extrusions.append(E1)

        tmp=extrusions[0]
        for n in range(1,len(extrusions)):
          tmp=tmp.fuse(extrusions[n])

        Part.show(tmp)
        FreeCAD.getDocument(FreeCAD.activeDocument().Name).ActiveObject.Label=plabel


# ******************************************************************************
# Extracts a board countour out of boardfile /plain/ section
# Currently not used / activated
def boardextract(brdfile, plabel, fusetol=0.0):
  extrusions = []
  edges = []
  
  board_outline_type = []
  board_outline_x1_or_Mx = []
  board_outline_y1_or_My = []
  board_outline_x2_or_r = []
  board_outline_y2 = []
  board_outline_curve = []
  board_outline_status = []
  
  LAYER_DIMENSION_NUMBER = "20"
  OUTLINE_TYPE_WIRE = 0
  OUTLINE_TYPE_ARC = 1
  OUTLINE_TYPE_CIRCLE = 2
  OUTLINE_STATUS_NEW  = 0
  OUTLINE_STATUS_USED = 1
  
  
  while len(edges)>0:
    edges.pop()
  while len(extrusions)>0:
    extrusions.pop()
  while len(board_outline_type)>0:
    board_outline_type.pop()
  while len(board_outline_x1_or_Mx)>0:
    board_outline_x1_or_Mx.pop()
  while len(board_outline_y1_or_My)>0:
    board_outline_y1_or_My.pop()
  while len(board_outline_x2_or_r)>0:
    board_outline_x2_or_r.pop()
  while len(board_outline_y2)>0:
    board_outline_y2.pop()
  while len(board_outline_curve)>0:
    board_outline_curve.pop()
    
  path1 = brdfile.getElementsByTagName("eagle")[0]
  path2 = path1.getElementsByTagName("drawing")[0]
  path3 = path2.getElementsByTagName("board")[0]
  path4 = path3.getElementsByTagName("plain")[0]
  wiresarcs = path4.getElementsByTagName("wire")
  circles = path4.getElementsByTagName("circle")
  
  for WA in wiresarcs:
    if WA.getAttribute("layer")==LAYER_DIMENSION_NUMBER:
      if WA.hasAttribute("curve"):
        board_outline_curve.append(float(WA.getAttribute("curve")))
        board_outline_type.append(OUTLINE_TYPE_ARC)
      else:
        board_outline_curve.append(float(0.0))
        board_outline_type.append(OUTLINE_TYPE_WIRE)
      board_outline_x1_or_Mx.append(float(WA.getAttribute("x1")))
      board_outline_y1_or_My.append(float(WA.getAttribute("y1")))
      board_outline_x2_or_r.append(float(WA.getAttribute("x2")))
      board_outline_y2.append(float(WA.getAttribute("y2")))
  
  for C in circles:
    if C.getAttribute("layer")==LAYER_DIMENSION_NUMBER:
      board_outline_curve.append(float(0.0))
      board_outline_type.append(OUTLINE_TYPE_CIRCLE)
      board_outline_x1_or_Mx.append(float(C.getAttribute("x")))
      board_outline_y1_or_My.append(float(C.getAttribute("y")))
      board_outline_x2_or_r.append(float(C.getAttribute("radius")))
      board_outline_y2.append(float(0.0))
  
  #testing
  #line   = Part.makeLine( beg, end )
  #arc    = Part.Edge( Part.Arc( Base.Vector(*beg), Base.Vector(*mid), Base.Vector(*end) ) )
  #helix  = Part.makeHelix( pitch, height, radius, degApexAngle )
  #circle = Part.makeCircle( radius, Base.Vector(*origin), Base.Vector(*direction) )
  #for OUT in board_outline_type:
  # todo: outline in lib nicht plain
  # todo outline circles
  for n in range(0,len(board_outline_type)):
    OUT = board_outline_type[n]
    if OUT==OUTLINE_TYPE_WIRE:
      FreeCAD.Console.PrintMessage("line\n")
      ln = Part.Line(Base.Vector(board_outline_x1_or_Mx[n],board_outline_y1_or_My[n],0.0),Base.Vector(board_outline_x2_or_r[n],board_outline_y2[n],0.0))
      wireseg = ln.toShape()
      #Part.show(wireseg)
      edges.append(wireseg)
    elif OUT==OUTLINE_TYPE_ARC:
      FreeCAD.Console.PrintMessage("arc\n")
      wx1 = board_outline_x1_or_Mx[n]
      wy1 = board_outline_y1_or_My[n]
      wx2 = board_outline_x2_or_r[n]
      wy2 = board_outline_y2[n]
      arccurve = board_outline_curve[n]
      ax1, ay1, a_middle_x, a_middle_y, ax2, ay2 = conv_arc_to_threepoint(wx1,wy1,wx2,wy2,arccurve)
      ar = Part.Arc(Base.Vector(ax1,ay1,0),Base.Vector(a_middle_x,a_middle_y,0),Base.Vector(ax2,ay2,0))
      wireseg = ar.toShape()
      #Part.show(wireseg)
      edges.append(wireseg)
  wi=Part.Wire(edges)
  #Part.show(W)
  if (wi.isClosed() == True):
    # face nur wenn true
    FACE = Part.Face(wi)
    EX = FACE.extrude(Base.Vector(0,0,10))
    Part.show(EX)
  else:
    Part.show(wi)
  #shape = Part.Shape(edges)
  #for v in shape.Vertexes: v.setTolerance(fusetol)
  #Part.show(shape)
  #W1 = Part.Wire(shape.Edges)
  #F1 = Part.Face(W1)
  #h =3.0
  #E1 = F1.extrude(Base.Vector(0,0,h))
  #Part.show(E1)
  #extrusions.append(E1)
  #tmp=extrusions[0]
  #for n in range(1,len(extrusions)):
  #tmp=tmp.fuse(extrusions[n])
  #Part.show(tmp)
  #FreeCAD.getDocument(FreeCAD.activeDocument().Name).ActiveObject.Label=plabel


# ******************************************************************************
# This function was the main functional block. It places the parts and calls 
# other functions for generationg the output
def placeparts(brdfile):
  #global pack_excludelist
  xmax=0.
  xmin=0.
  ymax=0.
  ymin=0.
  
  FreeCAD.newDocument(docname)
  
  # Parse the XML translation library file
  parselibfile()
  
  # Get the Eagle Board file elements
  path1 = brdfile.getElementsByTagName("eagle")[0]
  path2 = path1.getElementsByTagName("drawing")[0]
  path3 = path2.getElementsByTagName("board")[0]
  path4 = path3.getElementsByTagName("elements")[0]
  elements = path4.getElementsByTagName("element")
  
  for x in elements:
    name = x.getAttribute("name")
    package = x.getAttribute("package")
    posx = float(x.getAttribute("x"))
    posy = float(x.getAttribute("y"))
    
    # store min and max values for later use in PCB board-generation
    if (posx>xmax):
      xmax=posx
    elif (posx<xmin):
      xmin=posx

    if (posy>ymax):
      ymax=posy
    elif (posy<ymin):
      ymin=posy
    
    # Check if package of the part is in exclude list
    part_np=0
    for pex in pack_excludelist:
      # diagfile=__builtin__.open(diag_path+diag_filename, "a")
      # diagfile.writelines("exlist: "+pex+"\n")
      # diagfile.close()
      if pex==package:
        part_np=1
        break
      else:
        part_np=0

    if part_np==0:    
      angle_str = x.getAttribute("rot")
      if string.find(angle_str, "MR")==0:
        angle = float(string.lstrip(angle_str,"MR"))
        mirror = 1
      elif string.find(angle_str, "R")==0:
        angle = float(string.lstrip(angle_str,"R"))
        mirror = 0
      else:
        angle = 0.0
        mirror = 0
      
      try:
        idx = lib_pack.index(package)
        # Teil in Bibliothek gefunden also platzieren:
        #print lib_file[idx]
        if string.find(lib_file[idx], ".stp") != -1:
          Part.insert(libpath+lib_file[idx],docname)
          newinsert=1
        elif string.find(lib_file[idx], ".stl") != -1:
          Mesh.insert(libpath+lib_file[idx],docname)
          newinsert=1
        elif string.find(lib_file[idx], ".FCStd") != -1:
          insertfcstd(libpath+lib_file[idx],docname)
          newinsert=1
        else:
          newinsert=0
        
        if newinsert==1:
          # Setting the Part to its Zero-Offset position
          #FreeCAD.getDocument(docname).ActiveObject.Placement = FreeCAD.Placement(FreeCAD.Vector(float(lib_movx[idx]), float(lib_movy[idx]), float(lib_movz[idx])), FreeCAD.Rotation(eulerToQuaternion(float(lib_rotz[idx]), float(lib_roty[idx]), float(lib_rotx[idx]))))
          place = FreeCAD.Placement(FreeCAD.Vector(float(lib_movx[idx]), float(lib_movy[idx]), float(lib_movz[idx])), FreeCAD.Rotation(eulerToQuaternion(float(lib_rotz[idx]), float(lib_roty[idx]), float(lib_rotx[idx]))))
          FreeCAD.getDocument(docname).ActiveObject.Placement = place.multiply(FreeCAD.getDocument(docname).ActiveObject.Placement)
          
          if mirror==1:
            place = FreeCAD.Placement(FreeCAD.Vector(posx, posy, (-0.5*pcb_thickness)), FreeCAD.Rotation(eulerToQuaternion(360-angle,180,0)) )         
          else:
            place = FreeCAD.Placement(FreeCAD.Vector(posx, posy, (0.5*pcb_thickness)), FreeCAD.Rotation(eulerToQuaternion(angle,0,0)) )

          # Place the Part 
          FreeCAD.getDocument(docname).ActiveObject.Placement = place.multiply(FreeCAD.getDocument(docname).ActiveObject.Placement)
          # FreeCAD.getDocument(docname).ActiveObject.Label = name + "_" + FreeCAD.getDocument(docname).ActiveObject.Label
          FreeCAD.getDocument(docname).ActiveObject.Label = name + "_" + string.split(lib_file[idx],".")[0]
        
        if DIAGNOSTIC==True:
          diagfile=__builtin__.open(diag_path+diag_filename, "a")
          diagfile.writelines("Placed part: "+name+" with package "+package+"\n")
          diagfile.close()
      except ValueError:
        idx = 0
        if DIAGNOSTIC==True:
          diagfile=__builtin__.open(diag_path+diag_filename, "a")
          diagfile.writelines("Not found package: "+package+" for part "+name+"\n")
          diagfile.close()
        # Teil nicht in Bibliothek gefunden also generieren falls moeglich
        ##if partextract(....)==True:
        #partextract(dom, name + "_DUMMY", "RESC1608X55N", POINT_TOLERANCE)
        #if mirror==1:
        #  FreeCAD.getDocument(docname).ActiveObject.Placement = FreeCAD.Placement(FreeCAD.Vector(posx, posy, (-0.5*pcb_thickness)), FreeCAD.Rotation(eulerToQuaternion(360-angle, 180, 0)))
        #else:  
        #  FreeCAD.getDocument(docname).ActiveObject.Placement = FreeCAD.Placement(FreeCAD.Vector(posx, posy, (0.5*pcb_thickness)), FreeCAD.Rotation(eulerToQuaternion(angle, 0, 0)))
        if (show_notfound == PH_BALL):
          ballradius=0.25
          nf=Part.makeSphere(ballradius)
          if mirror==1:
            nf.translate(FreeCAD.Vector(posx, posy, ((-0.5*pcb_thickness)-ballradius)))
          else:
            nf.translate(FreeCAD.Vector(posx, posy, ((0.5*pcb_thickness)+ballradius)))
          Part.show(nf)
          FreeCADGui.getDocument(docname).ActiveObject.ShapeColor=(1.,0.,0.)
          FreeCAD.getDocument(docname).ActiveObject.Label = "NOTFOUND_" + name + "_" + package
    else:
      if DIAGNOSTIC==True:
        diagfile=__builtin__.open(diag_path+diag_filename, "a")
        diagfile.writelines("Not placed package: "+package+" for part "+name+" because of excludelist.\n")
        diagfile.close()

  # --- Place boardfile -------
  if (boardgentype == PLANEGEN):
    pcb1=Part.makeBox(xmax-xmin+bgen_border,ymax-ymin+bgen_border,pcb_thickness)
    pcb1.translate(FreeCAD.Vector( ((xmin+xmax)/2.)-((xmax-xmin+bgen_border)/2.), ((ymin+ymax)/2.)-((ymax-ymin+bgen_border)/2.), (-0.5*pcb_thickness)))
    Part.show(pcb1)
    FreeCADGui.getDocument(docname).ActiveObject.ShapeColor=(0.,(170./255.),0.)
    FreeCAD.getDocument(docname).ActiveObject.Label = "BOARD"
  elif (boardgentype == CADFILE):
    if string.find(boardCADfilePath, ".stp") != -1:
      Part.insert(boardCADfilePath,docname)
      newinsert=1
    elif string.find(boardCADfilePath, ".stl") != -1:
      Mesh.insert(boardCADfilePath,docname)
      newinsert=1
    else:
      newinsert=0

    if newinsert==1:
      #FreeCAD.getDocument(docname).ActiveObject.Placement = FreeCAD.Placement(FreeCAD.Vector(posx, posy, 0.5 * pcb_thickness), FreeCAD.Rotation(eulerToQuaternion(angle, 0, 0)))
      FreeCAD.getDocument(docname).ActiveObject.Label = "BOARD"
  
  # --- Place Signals ----------
  if (signalgen==True):
    showsignals(brdfile)
  
  # --- End of Placement -------
  FreeCAD.Gui.activeDocument().activeView().viewAxometric()
  FreeCAD.Gui.SendMsgToActiveView("ViewFit")


#placeparts(dom)
