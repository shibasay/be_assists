#!/usr/bin/env python
# fileencoding=utf-8
#author  : @shibasay
#date    : 2015/12/22

import sys
import re
import os
from optparse import OptionParser
ENV_HAS_PIL=None
try:
    from PIL import Image
    ENV_HAS_PIL = True
except ImportError:
    ENV_HAS_PIL = False

BE_VERSION=0.91

header=r"""#be! data file
#(int)x (int)y (int)z (int)r (int)g (int)b (int)a
#delimiters are space' ', tab'\t', comma',' and colon':'."""
header_new=r"""#be! data file
#At first,set field size
#(int)x (int)y (int)z
%d %d %d
#Next, set blocks
#(int)x (int)y (int)z (int)r (int)g (int)b (int)a
#delimiters are space' ', tab'\t', comma',' and colon':'."""

FILL_R=255
FILL_G=255
FILL_B=255
FILL_A=255

AXIS_X=0
AXIS_Y=1
AXIS_Z=2

axis_str_dic = {0: "X", 1: "Y", 2:"Z"}

import math
def rotate2D(x,y, degree):
    theta = math.radians(degree)
    cos = math.cos(theta)
    sin = math.sin(theta)
    # 浮動小数点->整数変換の誤差の影響を極力なくすため、
    # 90度の倍数の回転のときのみ sin(), cos() の結果を整数におとす
    # （理論上 1,-1,0 のどれかなので、演算前に整数にしても問題ない）
    if degree % 90 == 0:
        cos = int(cos)
        sin = int(sin)
    newx = cos * x - sin * y
    newy = sin * x + cos * y
    return newx, newy
def rotate3D(x,y,z, xdegree, ydegree, zdegree):
    newx,newy,newz = x,y,z
    newy, newz = rotate2D(newy,newz, xdegree)
    newx, newz = rotate2D(newx,newz, ydegree)
    newx, newy = rotate2D(newx,newy, zdegree)
    return newx, newy, newz
def ave(l): return reduce(lambda x,y:x+y, l, 0) / len(l)

class BED3D(dict):
    def __init__(self, seq=None, **kwargs):
        super(BED3D, self).__init__()
        self.xmax = None
        self.xmin = None
        self.ymax = None
        self.ymin = None
        self.zmax = None
        self.zmin = None
        self.xlistdic = None
        self.ylistdic = None
        self.zlistdic = None
        self.isPlane = False
        self.planeAxis = None # 0(x), 1(y), 2(z), None(not plane)
        self.planePos  = None
        self.planeWidth  = None
        self.planeHeight = None

    def printAll(self):
        outlist = []
        if BE_VERSION >= 0.9:
            self.updatePosMaxMin()
            outlist.append(header_new % (self.xmax, self.ymax, self.zmax))
        else:
            outlist.append(header)
        for v in self.values():
            outlist.append(v.getline())
        return "\n".join(outlist)

    def setIsPlane(self, axis, pos=0):
        self.isPlane = True
        self.planeAxis = axis
        self.planePos = pos
        self.updatePosMaxMin()
        if axis == AXIS_X:
            self.planeWidth  = self.ymax - self.ymin + 1
            self.planeHeight = self.zmax - self.zmin + 1
        elif axis == AXIS_Y:
            self.planeWidth  = self.xmax - self.xmin + 1
            self.planeHeight = self.zmax - self.zmin + 1
        elif axis == AXIS_Z:
            self.planeWidth  = self.xmax - self.xmin + 1
            self.planeHeight = self.ymax - self.ymin + 1

    def addBoxel(self, b):
        self[b.getPosTuple()] = b

    def getSurrounds(self, x,y,z):
        top  = self.get((x  ,y-1,z  ), None)
        bot  = self.get((x  ,y+1,z  ), None)
        lef  = self.get((x+1,y  ,z  ), None)
        rig  = self.get((x-1,y  ,z  ), None)
        fore = self.get((x  ,y  ,z+1), None)
        back = self.get((x  ,y  ,z-1), None)
        return top,bot,lef,rig,fore,back

    def updatePosMaxMin(self):
        xlist,ylist,zlist = zip(*self.keys()) # this is same to below
        #xlist = [x for x,y,z in self.keys()]
        #ylist = [y for x,y,z in self.keys()]
        #zlist = [z for x,y,z in self.keys()]
        self.xmax = max(xlist)
        self.xmin = min(xlist)
        self.ymax = max(ylist)
        self.ymin = min(ylist)
        self.zmax = max(zlist)
        self.zmin = min(zlist)

    def generateListXYis(self, xpos, ypos):
        return [self.get((xpos,ypos,z), None) for z in range(self.zmin, self.zmax+1)]
    def generateListYZis(self, ypos, zpos):
        return [self.get((x,ypos,zpos), None) for x in range(self.xmin, self.xmax+1)]
    def generateListZXis(self, zpos, xpos):
        return [self.get((xpos,y,zpos), None) for y in range(self.ymin, self.ymax+1)]

    def updateXYZlist(self):
        self.xlistdic = dict()
        self.ylistdic = dict()
        self.zlistdic = dict()
        for x in range(self.xmin, self.xmax+1):
            for y in range(self.ymin, self.ymax+1):
                self.zlistdic[x,y] = sorted([(e.z, e) for e in self.generateListXYis(x,y) if e])

        for z in range(self.zmin, self.zmax+1):
            for x in range(self.xmin, self.xmax+1):
                self.ylistdic[z,x] = sorted([(e.y, e) for e in self.generateListZXis(z,x) if e])

        for y in range(self.ymin, self.ymax+1):
            for z in range(self.zmin, self.zmax+1):
                self.xlistdic[y,z] = sorted([(e.x, e) for e in self.generateListYZis(y,z) if e])

    def getListXYis(self, x, y): return self.zlistdic[x,y]
    def getListYZis(self, y, z): return self.xlistdic[y,z]
    def getListZXis(self, z, x): return self.ylistdic[z,x]

    #============= out flag based delete decision ==========#
    def checkClosed_base_usingOut(self, value, level, smaller_list, larger_list):
        fand = lambda x,y: x and y
        smaller_adj_is_in = len(smaller_list) < level or reduce(fand, [not smaller_list[-i][1].outflag for i in range(1, level+1)], True)
        larger_adj_is_in = len(larger_list) < level or reduce(fand, [not larger_list[i-1][1].outflag for i in range(1, level+1)], True)
        return (smaller_adj_is_in and larger_adj_is_in)

    def checkClosedX_usingOut(self, x,y,z,  level):
        xlist = self.getListYZis(y,z)
        xlist_roi = [(ex, e) for (ex, e) in xlist if x-level <= ex <= x+level]
        idx = xlist_roi.index((x, self[x,y,z])) # find self
        xlist_smaller = xlist_roi[:idx]
        xlist_larger  = xlist_roi[idx+1:]
        return self.checkClosed_base_usingOut(x, level, xlist_smaller, xlist_larger)

    def checkClosedY_usingOut(self, x,y,z, level):
        ylist = self.getListZXis(z,x)
        ylist_roi = [(ey, e) for (ey, e) in ylist if y-level <= ey <= y+level]
        idx = ylist_roi.index((y, self[x,y,z])) # find self
        ylist_smaller = ylist_roi[:idx]
        ylist_larger  = ylist_roi[idx+1:]
        return self.checkClosed_base_usingOut(y, level, ylist_smaller, ylist_larger)

    def checkClosedZ_usingOut(self, x,y,z, level):
        zlist = self.getListXYis(x,y)
        zlist_roi = [(ez, e) for (ez, e) in zlist if z-level <= ez <= z+level]
        idx = zlist_roi.index((z, self[x,y,z])) # find self
        zlist_smaller = zlist_roi[:idx]
        zlist_larger  = zlist_roi[idx+1:]
        return self.checkClosed_base_usingOut(z, level, zlist_smaller, zlist_larger)

    def checkClosedPos_usingOut(self, x,y,z, level):
        allclosed = (self.checkClosedX_usingOut(x,y,z, level) and
                     self.checkClosedY_usingOut(x,y,z, level) and
                     self.checkClosedZ_usingOut(x,y,z, level))
        return allclosed

    def initOutFlag(self):
        for v in self.values():
            v.outflag = False
            v.doneflag = False

    def cleanOuts(self):
        newmodel = BED3D()
        for k, b in self.items():
            if not b.outflag:
                newmodel.addBoxel(b)
        return newmodel

    def getFillingOut(self, x,y,z):
        xin = self.xmin <= x <= self.xmax+1
        yin = self.ymin <= y <= self.ymax+1
        zin = self.zmin <= z <= self.zmax+1
        if xin and yin and zin:
            b = self.get((x,y,z), None)
            if b is None:
                b = Boxel(x, y, z, FILL_R, FILL_G, FILL_B, FILL_A)
                b.outflag = True
                self.addBoxel(b)
            return b
        else:
            return None

    def getSurroundsFillingOut(self, x,y,z):
        top  = self.getFillingOut(x  , y+1, z  )
        bot  = self.getFillingOut(x  , y-1, z  )
        rig  = self.getFillingOut(x+1, y  , z  )
        lef  = self.getFillingOut(x-1, y  , z  )
        back = self.getFillingOut(x  , y  , z+1)
        fore = self.getFillingOut(x  , y  , z-1)
        return top,bot,lef,rig,fore,back

    def getSurroundingUndoneOuts(self,b):
        b.doneflag = True
        surrounds = self.getSurroundsFillingOut(b.x, b.y, b.z)
        return [b for b in surrounds if b and b.outflag and not b.doneflag]

    def setOutBoxels(self):
        self.initOutFlag()
        self.updatePosMaxMin()

        newmodel = BED3D()
        for k,b in self.items(): newmodel.addBoxel(b)

        # make base position which is absolutely OUT
        n = Boxel(self.xmin-1, self.ymin-1, self.zmin-1,FILL_R, FILL_G, FILL_B, FILL_A)
        n.outflag = True
        newmodel.addBoxel(n)
        newmodel.updatePosMaxMin()

        undones = newmodel.getSurroundingUndoneOuts(n)
        while undones:
            newundones = set()
            for b in undones:
                for e in newmodel.getSurroundingUndoneOuts(b):
                    newundones.add(e)
            undones = newundones
        newmodel.updatePosMaxMin()
        return newmodel

    def delClosed_remainingOut(self, level):
        outfilled = self.setOutBoxels()
        outfilled.updateXYZlist() # prepare checkClosed
        newmodel = BED3D()
        for (x,y,z), v in outfilled.items():
            if not outfilled.checkClosedPos_usingOut(x,y,z, level):
                newmodel.addBoxel(v)
        return newmodel

    def delClosed(self, level):
        return self.delClosed_remainingOut(level).cleanOuts()

    #==== Fill closed regions with boxels ====#
    def fillClosed_remainingOut(self):
        outfilled = self.setOutBoxels()
        newmodel = BED3D()
        for x in range(outfilled.xmin, outfilled.xmax+1):
            for y in range(outfilled.ymin, outfilled.ymax+1):
                for z in range(outfilled.zmin, outfilled.zmax+1):
                    #print "check (%d, %d, %d)" % (x,y,z),
                    v = outfilled.get((x,y,z), None)
                    if v: # exist original (or out) boxel
                        newmodel.addBoxel(v)
                    else: # do not exist boxel, so fill there
                        surrounds_rgb = [(e.r, e.g, e.b) for e in self.getSurrounds(x,y,z) if e and not e.outflag]
                        if surrounds_rgb:
                            rlist,glist,blist = zip(*surrounds_rgb)
                            newmodel.addBoxel(Boxel(x,y,z, ave(rlist), ave(glist), ave(blist), FILL_A))
                        else:
                            newmodel.addBoxel(Boxel(x,y,z, FILL_R, FILL_G, FILL_B, FILL_A))
        return newmodel

    def fillClosed(self):
        return self.fillClosed_remainingOut().cleanOuts()

    #==== Mirror model generation ====#
    def makeMirrorModel(self, mode):
        # entry
        if mode==0:   # center mirror (use left)
            newmodel, oddflag = self.delHalf(0)
            return newmodel.makeMirror(0, oddflag)
        elif mode==1: # center mirror (use right)
            newmodel, oddflag = self.delHalf(1)
            return newmodel.makeMirror(1, oddflag)
        else:         # whole mirror (make other one)
            return self.makeMirror(0, 0)

    def delHalf(self, lr):
        self.updatePosMaxMin()
        width = self.xmax - self.xmin + 1
        offset = width & 1
        center = width//2 + offset - 1 + self.xmin
        newmodel = BED3D()
        #print "@delHalf: width = %d, offset = %d, center = %d, xmin = %d, xmax = %d" % (width, offset, center, self.xmin, self.xmax)
        if lr == 0: # delete right
            for (x,y,z), v in self.items():
                if x <= center:
                    newmodel.addBoxel(v)
        else:
            for (x,y,z), v in self.items():
                if x > center-offset:
                    newmodel.addBoxel(v)
        return newmodel, offset

    def makeMirror(self, lr, oddflag):
        self.updatePosMaxMin()
        center = self.xmax if lr==0 else self.xmin
        width = self.xmax - self.xmin + 1
        offset = 1 - oddflag
        newmodel = BED3D()
        #print "@makeMirror: width = %d, oddflag = %d, xmin = %d, xmax = %d" % (width, oddflag, self.xmin, self.xmax)
        for (x,y,z), v in self.items():
            newmodel.addBoxel(v)
            diff = x-center
            newx = center-diff+offset if lr==0 else center-diff-offset
            if newx != x:
                newmodel.addBoxel(v.genModPosBoxel(newx,y,z))
        return newmodel

    #==== Scaled model generation ====#
    def makeScaledModel(self, scale):
        newmodel = self.positionScaledBoxels(scale)
        newmodel = newmodel.genScaledBoxels(scale)
        newmodel = newmodel.smoothing(scale)
        return newmodel

    def positionScaledBoxels(self, scale):
        newmodel = BED3D()
        for (x,y,z), b in self.items():
            scalex = x * scale
            scaley = y * scale
            scalez = z * scale
            newmodel.addBoxel(b.genModPosBoxel(scalex, scaley, scalez))
        return newmodel

    def genScaledBoxels(self, scale):
        newmodel = BED3D()
        for (x,y,z), b in self.items():
            for newz in range(z, z+scale):
                for newy in range(y, y+scale):
                    for newx in range(x, x+scale):
                        newmodel.addBoxel(b.genModPosBoxel(newx, newy, newz))
        return newmodel

    def smoothing(self, scale):
        outfilled = self.fillClosed_remainingOut()
        originals = [((x,y,z), b) for (x,y,z), b in self.items() if x%scale==0 and y%scale==0 and z%scale==0]
        newmodel = BED3D()
        for (x,y,z),b in originals:
            scales = []
            for zz in range(z, z+scale):
                for yy in range(y, y+scale):
                    for xx in range(x, x+scale):
                        scales.append(self[(xx,yy,zz)])
            surrounds = []
            for zz in range(z, z+scale):
                for yy in range(y, y+scale):
                    surrounds.append(outfilled[(x-1,yy,zz)])
                    surrounds.append(outfilled[(x+scale,yy,zz)])
            for yy in range(y, y+scale):
                for xx in range(x, x+scale):
                    surrounds.append(outfilled[(xx,yy,z-1)])
                    surrounds.append(outfilled[(xx,yy,z+scale)])
            for xx in range(x, x+scale):
                for zz in range(z, z+scale):
                    surrounds.append(outfilled[(xx,y-1,zz)])
                    surrounds.append(outfilled[(xx,y+scale,zz)])
            surrounds_in = [e for e in surrounds if not e.outflag]

            for s in scales:
                if not set(outfilled.getSurrounds(s.x, s.y, s.z)) & set(surrounds_in):
                    # boxels surfacing out is conditionally deleted
                    surface_outs = [e for e in outfilled.getSurrounds(s.x, s.y, s.z) if e.outflag]
                    if len(surface_outs) < 3:
                        newmodel.addBoxel(s)
                else:
                    # boxels surfacing inside always remain
                    newmodel.addBoxel(s)
        return newmodel

    #==== moved model generation ====#
    def makeMovedModel(self, mvx, mvy, mvz):
        newmodel = BED3D()
        for (x,y,z), b in self.items():
            newmodel.addBoxel(b.genModPosBoxel(x+mvx, y+mvy, z+mvz))
        return newmodel

    #==== rotated model generation ====#
    def makeRotatedModel(self, xdegree, ydegree, zdegree):
        newmodel = BED3D()
        for (x,y,z), b in self.items():
            newx, newy, newz = [int(v) for v in rotate3D(x,y,z, xdegree, ydegree, zdegree)]
            newmodel.addBoxel(b.genModPosBoxel(newx, newy, newz))
        return newmodel

    #==== slim/fat model generation ====#
    def makeSplitModels(self, axis, pos):
        # axis = 0:x, 1:y, 2:z
        newmodel1 = BED3D()
        newmodel2 = BED3D()
        newmodel3 = BED3D()
        for (x,y,z), b in self.items():
            if axis == AXIS_X:
                if x <  pos: newmodel1.addBoxel(b.genModPosBoxel(x, y, z))
                if x == pos: newmodel2.addBoxel(b.genModPosBoxel(x, y, z))
                if x >  pos: newmodel3.addBoxel(b.genModPosBoxel(x, y, z))
            elif axis == AXIS_Y:
                if y <  pos: newmodel1.addBoxel(b.genModPosBoxel(x, y, z))
                if y == pos: newmodel2.addBoxel(b.genModPosBoxel(x, y, z))
                if y >  pos: newmodel3.addBoxel(b.genModPosBoxel(x, y, z))
            elif axis == AXIS_Z:
                if z <  pos: newmodel1.addBoxel(b.genModPosBoxel(x, y, z))
                if z == pos: newmodel2.addBoxel(b.genModPosBoxel(x, y, z))
                if z >  pos: newmodel3.addBoxel(b.genModPosBoxel(x, y, z))
            else:
                raise "unsupported value for axis: %d" % axis
        newmodel2.setIsPlane(axis, pos)
        return newmodel1, newmodel2, newmodel3

    def makeJoinedModel(self, another):
        # just generate a new model by joining two models without any position modification
        newmodel = BED3D()
        for (x,y,z), b in self.items():
            newmodel.addBoxel(b.genModPosBoxel(x,y,z))
        for (x,y,z), b in another.items():
            newmodel.addBoxel(b.genModPosBoxel(x,y,z))
        return newmodel

    def makeSlimModel(self, axis, pos):
        # delete 1 plane specified by using axis and pos
        m1, m2, m3 = self.makeSplitModels(axis, pos)
        if   axis == AXIS_X: (mvx, mvy, mvz) = (-1,  0,  0)
        elif axis == AXIS_Y: (mvx, mvy, mvz) = ( 0, -1,  0)
        elif axis == AXIS_Z: (mvx, mvy, mvz) = ( 0,  0, -1)
        else: raise "unsupported value for axis: %d" % axis
        return m1.makeJoinedModel(m3.makeMovedModel(mvx,mvy,mvz))

    def makeFatModel(self, axis, pos):
        # duplicate 1 plane specified by using axis and pos
        m1, m2, m3 = self.makeSplitModels(axis, pos)
        if   axis == AXIS_X: (mvx, mvy, mvz) = ( 1,  0,  0)
        elif axis == AXIS_Y: (mvx, mvy, mvz) = ( 0,  1,  0)
        elif axis == AXIS_Z: (mvx, mvy, mvz) = ( 0,  0,  1)
        else: raise "unsupported value for axis: %d" % axis
        m2a = m2.makeMovedModel(mvx,mvy,mvz)
        m3a = m3.makeMovedModel(mvx,mvy,mvz)
        return m1.makeJoinedModel(m2).makeJoinedModel(m2a).makeJoinedModel(m3a)

    #==== projected model generation ====#
    def makeProjectedModels(self, axis, minmaxflag):
        self.updateXYZlist()
        newmodel = BED3D()
        if axis == AXIS_X: # project to YZ plane (reduce X axis)
            for y in range(self.ymin, self.ymax+1):
                for z in range(self.zmin, self.zmax+1):
                    l = self.getListYZis(y, z)
                    if l:
                        pos, b = min(l) if minmaxflag == 0 else max(l)
                        newmodel.addBoxel(b.genModPosBoxel(0,y,z))
        elif axis == AXIS_Y: # project to XZ plane (reduce Y axis)
            for z in range(self.zmin, self.zmax+1):
                for x in range(self.xmin, self.xmax+1):
                    l = self.getListZXis(z, x)
                    if l:
                        pos, b = min(l) if minmaxflag == 0 else max(l)
                        newmodel.addBoxel(b.genModPosBoxel(x,0,z))
        elif axis == AXIS_Z: # project to XY plane (reduce Z axis)
            for x in range(self.xmin, self.xmax+1):
                for y in range(self.ymin, self.ymax+1):
                    l = self.getListXYis(x, y)
                    if l:
                        pos, b = min(l) if minmaxflag == 0 else max(l)
                        newmodel.addBoxel(b.genModPosBoxel(x,y,0))
        else:
            raise "unexptected axis %d is provided..." % axis
        newmodel.setIsPlane(axis)
        return newmodel

    def generatePNG(self, outnamebase, voxelsize=1, path="."):
        assert self.isPlane, "generatePNG should be called with plane model!"
        img = Image.new("RGBA", (self.planeWidth*voxelsize, self.planeHeight*voxelsize))
        for b in self.values():
            x,y, r,g,b,a = b.getPlaneData(self.planeAxis, self.xmin, self.ymin, self.zmin)
            for i in range(voxelsize):
                for j in range(voxelsize):
                    img.putpixel((x*voxelsize+i,y*voxelsize+j), (r,g,b,a))
        axis_str = None

        outname = "%s_%s.png" % (outnamebase, axis_str_dic[self.planeAxis])
        outpath = os.path.join(path, outname)
        img.save(outpath)
        return outpath

    def getPreviewXY_path(self, basename, path):
        self.updatePosMaxMin()
        newmodel = self.makeProjectedModels(AXIS_Z, 0)
        return newmodel.generatePNG(basename, 4, path)

    def getPreviewYZ_path(self, basename, path):
        self.updatePosMaxMin()
        newmodel = self.makeProjectedModels(AXIS_X, 0)
        return newmodel.generatePNG(basename, 4, path)

    def getPreviewZX_path(self, basename, path):
        self.updatePosMaxMin()
        newmodel = self.makeProjectedModels(AXIS_Y, 0)
        return newmodel.generatePNG(basename, 4, path)

    def makeFlippedModel(self, axis):
        self.updatePosMaxMin()
        newmodel = BED3D()
        for (x,y,z), b in self.items():
            if   axis == AXIS_X: newx, newy, newz = self.xmax - x, y, z
            elif axis == AXIS_Y: newx, newy, newz = x, self.ymax - y, z
            elif axis == AXIS_Z: newx, newy, newz = x, y, self.zmax - z
            newmodel.addBoxel(b.genModPosBoxel(newx, newy, newz))
        return newmodel

class Boxel(object):
    def __init__(self, x,y,z, r,g,b, alpha):
        self.x = int(x)
        self.y = int(y)
        self.z = int(z)
        self.r = int(r)
        self.g = int(g)
        self.b = int(b)
        self.alpha = int(alpha)
        self.outflag = False
        self.doneflag = False

    def getline(self):
        return "%d %d %d %d %d %d %d" % (self.x,self.y,self.z, self.r, self.g, self.b, self.alpha)

    def getPosTuple(self):
        return (self.x, self.y, self.z)

    def __repr__(self):
        return self.getline()

    def genModPosBoxel(self,x,y,z):
        newbed = Boxel(x, y, z, self.r, self.g, self.b, self.alpha)
        return newbed

    def getPlaneData(self, axis, xmin, ymin, zmin):
        if   axis == AXIS_X: return self.y-ymin, self.z-zmin, self.r, self.g, self.b, self.alpha
        elif axis == AXIS_Y: return self.x-xmin, self.z-zmin, self.r, self.g, self.b, self.alpha
        elif axis == AXIS_Z: return self.x-xmin, self.y-ymin, self.r, self.g, self.b, self.alpha

def bed_read(bedfilename):
    bed3D = BED3D()
    with open(bedfilename) as f:
        for line in f:
            if line[0] != '#':
                x,y,z, r,g,b, alpha = re.split(r"[ ,;\t]", line.strip())
                b = Boxel(x,y,z, r,g,b, alpha)
                bed3D[b.getPosTuple()] = b
    return bed3D

def getopt():
    version = '%prog 0.1'
    # usageの %prog はOptionParserによってos.path.basename(sys.argv[0])に置換えられる
    parser = OptionParser(usage=None, version=version)
    parser.add_option("-i", "--input",
                      dest="inputfile",
                      help="specify input file name",
                      metavar="INPUT")
    parser.add_option("-o", "--output",
                      dest="outputfile",
                      help="specify output file name",
                      metavar="OUTPUT")
    parser.add_option("-m", "--mode",
                      dest="mode",
                      type="int",
                      default=9999,
                      help="specify processing mode: 0(hollow), 1(fill), 2(mirror), 3(scale), 4(move), 5(rotate), 6(slim), 7(fat), 9999(info)",
                      metavar="MODE")
    # options for hollowing
    parser.add_option("-l", "--level",
                      dest="level",
                      type="int",
                      default=1,
                      help="specify level for hollowing mode: should be greater than or equal to 0",
                      metavar="LEVEL")
    # options for mirror generation
    parser.add_option("-r", "--mirrormode",
                      dest="mirrormode",
                      type="int",
                      default=0,
                      help="specify mirrormode: 0(use left half), 1(use right half), 2(make other one)",
                      metavar="LEVEL")
    # options for scaling
    parser.add_option("-s", "--scale",
                      dest="scale",
                      type="int",
                      default=1,
                      help="specify scaling factor integer",
                      metavar="SCALE")
    # options for move
    parser.add_option("-x", "--mvx",
                      dest="mvx",
                      type="int",
                      default=0,
                      help="specify move vector x",
                      metavar="MVX")
    parser.add_option("-y", "--mvy",
                      dest="mvy",
                      type="int",
                      default=0,
                      help="specify move vector y",
                      metavar="MVY")
    parser.add_option("-z", "--mvz",
                      dest="mvz",
                      type="int",
                      default=0,
                      help="specify move vector z",
                      metavar="MVZ")
    # options for rotation
    parser.add_option("--xdegree",
                      dest="xdegree",
                      type="int",
                      default=0,
                      help="specify rotation degree on axis x",
                      metavar="XDEGREE")
    parser.add_option("--ydegree",
                      dest="ydegree",
                      type="int",
                      default=0,
                      help="specify rotation degree on axis y",
                      metavar="YDEGREE")
    parser.add_option("--zdegree",
                      dest="zdegree",
                      type="int",
                      default=0,
                      help="specify rotation degree on axis z",
                      metavar="ZDEGREE")
    # options for slim/fat generation and projected model generation
    parser.add_option("-a", "--axis",
                      dest="axis",
                      type="int",
                      default=0,
                      help="specify axis for slim/fat plane position: 0(x), 1(y), 2(z)",
                      metavar="AXIS")
    parser.add_option("-p", "--p",
                      dest="pos",
                      type="int",
                      default=0,
                      help="specify plane position",
                      metavar="POS")
    # options for projected model generation
    parser.add_option("--minmaxflag",
                      dest="minmaxflag",
                      type="int",
                      default=0,
                      help="specify projection direction min or max",
                      metavar="MINMAX")
    (options, args) = parser.parse_args() # 引数パーズ

    if options.inputfile is None:
        parser.print_help()
        exit()
    else:
        return options


if __name__ == "__main__":
    options = getopt()

    bed3D = bed_read(options.inputfile)
    bed3D.updatePosMaxMin()
    sys.stderr.write("input model size info:\n")
    sys.stderr.write("\tmin  (x,y,z) = ( %3d, %3d, %3d)\n" % (bed3D.xmin, bed3D.ymin, bed3D.zmin))
    sys.stderr.write("\tmax  (x,y,z) = ( %3d, %3d, %3d)\n" % (bed3D.xmax, bed3D.ymax, bed3D.zmax))
    sys.stderr.write("\tdiff (x,y,z) = ( %3d, %3d, %3d)\n" % (bed3D.xmax-bed3D.xmin, bed3D.ymax-bed3D.ymin, bed3D.zmax-bed3D.zmin))

    if options.mode == 9999: exit()

    exedict = {
        0: lambda options: bed3D.delClosed(options.level),
        1: lambda options: bed3D.fillClosed(),
        2: lambda options: bed3D.makeMirrorModel(options.mirrormode),
        3: lambda options: bed3D.makeScaledModel(options.scale),
        4: lambda options: bed3D.makeMovedModel(options.mvx, options.mvy, options.mvz),
        5: lambda options: bed3D.makeRotatedModel(options.xdegree, options.ydegree, options.zdegree),
        6: lambda options: bed3D.makeSlimModel(options.axis, options.pos),
        7: lambda options: bed3D.makeFatModel(options.axis, options.pos),
        8: lambda options: bed3D.makeProjectedModels(options.axis, options.minmaxflag),
        }
    newmodel = exedict[options.mode](options)
    outstr = newmodel.printAll()

    if options.outputfile:
        with open(options.outputfile, "wt") as of:
            of.write(outstr)
    else:
        sys.stdout.write(outstr)

    if ENV_HAS_PIL and newmodel.isPlane:
        outnamebase = options.outputfile[:options.outputfile.rfind(".")]
        newmodel.generatePNG(outnamebase)

