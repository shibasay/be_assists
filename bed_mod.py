#!/usr/bin/env python
#author  : @shibasay 
#date    : 2014/5/5

from optparse import OptionParser

header=r"""#be! data file
#(int)x (int)y (int)z (int)r (int)g (int)b (int)a
#delimiters are space' ', tab'\t', camma',' and colon':'."""

FILL_R=255
FILL_G=255
FILL_B=255
FILL_A=255

class BED3D(dict):
    def printAll(self):
        outlist = []
        outlist.append(header)
        for v in self.values():
            outlist.append(v.getline())
        return "\n".join(outlist)

#    def getListXis(self, xpos):
#        return [v for (x,y,z), v in self.items() if x==xpos]
#    def getListYis(self, ypos):
#        return [v for (x,y,z), v in self.items() if y==ypos]
#    def getListZis(self, zpos):
#        return [v for (x,y,z), v in self.items() if z==zpos]

    def getSurrounds(self, x,y,z): 
        top = self.get((x,y-1,z), None)
        bot = self.get((x,y+1,z), None)
        lef = self.get((x+1,y,z), None)
        rig = self.get((x-1,y,z), None)
        fore = self.get((x,y,z+1), None)
        back = self.get((x,y,z-1), None)
        return top,bot,lef,rig,fore,back

    def updatePosMaxMin(self):
        xlist = [x for x,y,z in self.keys()]
        ylist = [y for x,y,z in self.keys()]
        zlist = [z for x,y,z in self.keys()]
        self.xmax = max(xlist)
        self.xmin = min(xlist)
        self.ymax = max(ylist)
        self.ymin = min(ylist)
        self.zmax = max(zlist)
        self.zmin = min(zlist)

    def generateListXYis(self, xpos, ypos):
        return [v for (x,y,z), v in self.items() if x==xpos and y==ypos]
    def generateListYZis(self, ypos, zpos):
        return [v for (x,y,z), v in self.items() if y==ypos and z==zpos]
    def generateListZXis(self, zpos, xpos):
        return [v for (x,y,z), v in self.items() if z==zpos and x==xpos]

    def updateXYZlist(self):
        self.xlistdic = dict()
        self.ylistdic = dict()
        self.zlistdic = dict()
        for x in range(self.xmin, self.xmax+1):
            for y in range(self.ymin, self.ymax+1):
                self.zlistdic[x,y]  = sorted([(e.z, e) for e in self.generateListXYis(x,y)])

        for z in range(self.zmin, self.zmax+1):
            for x in range(self.xmin, self.xmax+1):
                self.ylistdic[z,x]  = sorted([(e.y, e) for e in self.generateListZXis(z,x)])

        for y in range(self.ymin, self.ymax+1):
            for z in range(self.zmin, self.zmax+1):
                self.xlistdic[y,z]  = sorted([(e.x, e) for e in self.generateListYZis(y,z)])

    def getListXYis(self, x, y):
        return self.zlistdic[x,y]
    def getListYZis(self, y, z):
        return self.xlistdic[y,z]
    def getListZXis(self, z, x):
        return self.ylistdic[z,x]

    #============= out flag based delete decision ==========#
    def checkClosed_base_usingOut(self, value, level, smaller_list, larger_list):
        def fand(x,y): return (x and y)
        closed_weak = len(smaller_list) >= level and len(larger_list) >= level
        smaller_adj_is_in = True
        larger_adj_is_in = True
        if closed_weak:
            smaller_adj_is_in = reduce(fand, [(smaller_list[-i][0] != value-i or smaller_list[-i][1].outflag != True) for i in range(1, level+1)], True)
            larger_adj_is_in = reduce(fand, [(larger_list[i-1][0] != value+i or larger_list[i-1][1].outflag != True) for i in range(1, level+1)], True)
        closed = closed_weak and smaller_adj_is_in and larger_adj_is_in
        return closed
    def checkClosedX_usingOut(self,x,y,z,  level):
        xlist = self.getListYZis(y,z)
        xlist_smaller = [(ex, e) for (ex, e) in xlist if ex < x]
        xlist_larger  = [(ex, e) for (ex, e) in xlist if ex > x]
        return self.checkClosed_base_usingOut(x, level, xlist_smaller, xlist_larger)
    def checkClosedY_usingOut(self,x,y,z, level):
        ylist = self.getListZXis(z,x)
        ylist_smaller = [(ey, e) for (ey, e) in ylist if ey < y]
        ylist_larger  = [(ey, e) for (ey, e) in ylist if ey > y]
        return self.checkClosed_base_usingOut(y, level, ylist_smaller, ylist_larger)
    def checkClosedZ_usingOut(self,x,y,z, level):
        zlist = self.getListXYis(x,y)
        zlist_smaller = [(ez, e) for (ez, e) in zlist if ez < z]
        zlist_larger  = [(ez, e) for (ez, e) in zlist if ez > z]
        return self.checkClosed_base_usingOut(z, level, zlist_smaller, zlist_larger)
    def checkClosedPos_usingOut(self,x,y,z, level):
        allclosed = (self.checkClosedX_usingOut(x,y,z, level) and
                     self.checkClosedY_usingOut(x,y,z, level) and
                     self.checkClosedZ_usingOut(x,y,z, level))
        return allclosed

    def initOutFlag(self):
        for v in self.values():
            v.outflag = None
            v.doneflag = False

    def cleanOuts(self):
        newmodel = BED3D()
        for k, b in self.items():
            if b.outflag != True:
                newmodel[k] = b
        return newmodel

    def getFillingOut(self,x,y,z):
        if (self.xmin <= x <= self.xmax+1 and 
            self.ymin <= y <= self.ymax+1 and
            self.zmin <= z <= self.zmax+1):
            b = self.get((x,y,z), None)
            if b == None:
                b = Boxel(x, y, z, FILL_R, FILL_G, FILL_B, FILL_A)
                b.outflag = True
                self[(x,y,z)] = b
            else:
                if b.outflag == None:
                    b.outflag = False
            return b
        else:
            return None
    def getSurroundsFillingOut(self, x,y,z):
        top  = self.getFillingOut(x,y+1,z)
        bot  = self.getFillingOut(x,y-1,z)
        rig  = self.getFillingOut(x+1,y,z)
        lef  = self.getFillingOut(x-1,y,z)
        back = self.getFillingOut(x,y,z+1)
        fore = self.getFillingOut(x,y,z-1)
        return top,bot,lef,rig,fore,back
    def getSurroundingUndoneOuts(self,b):
        b.doneflag = True
        surrounds = self.getSurroundsFillingOut(b.x, b.y, b.z)
        return [b for b in surrounds if b and b.outflag == True and b.doneflag != True]
    def setOutBoxels(self):
        self.initOutFlag()
        self.updatePosMaxMin()

        newdata = BED3D()
        for k,b in self.items(): newdata[k] = b

        # make base position which is absolutely OUT
        n = Boxel(self.xmin-1, self.ymin-1, self.zmin-1,FILL_R, FILL_G, FILL_B, FILL_A)
        n.outflag = True
        newdata[n.getPosTuple()] = n
        newdata.updatePosMaxMin()

        undones = newdata.getSurroundingUndoneOuts(n)
        while undones:
            newundones = set()
            for b in undones:
                for e in newdata.getSurroundingUndoneOuts(b):
                    newundones.add(e)
            undones = newundones
        newdata.updatePosMaxMin()
        return newdata

    def delClosed_remainingOut(self, level):
        outfilled = self.setOutBoxels()
        outfilled.updateXYZlist()
        newdata = BED3D()
        for (x,y,z), v in outfilled.items():
            if not outfilled.checkClosedPos_usingOut(x,y,z, level):
                newdata[(x,y,z)] = v
        return newdata
    def delClosed(self, level):
        return self.delClosed_remainingOut(level).cleanOuts()

    def fillClosed_remainingOut(self):
        outfilled = self.setOutBoxels()
        newdata = BED3D()
        for x in range(outfilled.xmin, outfilled.xmax+1):
            for y in range(outfilled.ymin, outfilled.ymax+1):
                for z in range(outfilled.zmin, outfilled.zmax+1):
                    #print "check (%d, %d, %d)" % (x,y,z), 
                    v = outfilled.get((x,y,z), None)
                    if v:
                        newdata[(x,y,z)] = v
                        #print "exist original data" 
                    else:
                        #fillboxel = "%d,%d,%d %d,%d,%d %d" % (x,y,z, FILL_R, FILL_G, FILL_B, FILL_A)
                        n = Boxel(x,y,z, FILL_R, FILL_G, FILL_B, FILL_A)
                        newdata[n.getPosTuple()] = n
                        #print "fill!" 
        return newdata
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
        newdata = BED3D()
        #print "@delHalf: width = %d, offset = %d, center = %d, xmin = %d, xmax = %d" % (width, offset, center, self.xmin, self.xmax)
        if lr == 0: # delete right
            for (x,y,z), v in self.items():
                if x <= center:
                    newdata[(x,y,z)] = v
        else:
            for (x,y,z), v in self.items():
                if x > center-offset:
                    newdata[(x,y,z)] = v
        return newdata, offset

    def makeMirror(self, lr, oddflag):
        self.updatePosMaxMin()
        center = self.xmax if lr==0 else self.xmin
        width = self.xmax - self.xmin + 1
        offset = 1 - oddflag
        newdata = BED3D()
        #print "@makeMirror: width = %d, oddflag = %d, xmin = %d, xmax = %d" % (width, oddflag, self.xmin, self.xmax)
        for (x,y,z), v in self.items():
            newdata[(x,y,z)] = v
            diff = x-center
            newx = center-diff+offset if lr==0 else center-diff-offset
            if newx != x:
                newdata[(newx,y,z)] = v.genModPos(newx,y,z)
        return newdata

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
            newmodel[(scalex, scaley, scalez)] = b.genModPos(scalex, scaley, scalez)
        return newmodel
    def genScaledBoxels(self, scale):
        newmodel = BED3D()
        for (x,y,z), b in self.items():
            for newz in  range(z, z+scale):
                for newy in  range(y, y+scale):
                    for newx in  range(x, x+scale):
                        newmodel[(newx, newy, newz)] = b.genModPos(newx, newy, newz)
        return newmodel
    def smoothing(self, scale):
        outfilled = self.fillClosed_remainingOut()
        newmodel = BED3D()
        originals = [((x,y,z), b) for (x,y,z), b in self.items() if x%scale==0 and y%scale==0 and z%scale==0]
        for (x,y,z),b in originals:
            #print "orignal (%d, %d, %d)" % (x,y,z)
            scales = []
            for zz in  range(z, z+scale):
                for yy in  range(y, y+scale):
                    for xx in  range(x, x+scale):
                        scales.append(self[(xx,yy,zz)])
                        #print "append scales = (x,y,z) = (%d, %d, %d)" % (xx, yy, zz)
            surrounds = []
            for zz in  range(z, z+scale):
                for yy in  range(y, y+scale):
                    surrounds.append(outfilled[(x-1,yy,zz)])
                    surrounds.append(outfilled[(x+scale,yy,zz)])
            for yy in  range(y, y+scale):
                for xx in  range(x, x+scale):
                    surrounds.append(outfilled[(xx,yy,z-1)])
                    surrounds.append(outfilled[(xx,yy,z+scale)])
            for xx in  range(x, x+scale):
                for zz in  range(z, z+scale):
                    surrounds.append(outfilled[(xx,y-1,zz)])
                    surrounds.append(outfilled[(xx,y+scale,zz)])
            surrounds_in = [e for e in surrounds if e.outflag != True]
            #print "check scales"
            for s in scales:
                #print "\t(x,y,z) = (%d, %d, %d)" % (s.x, s.y, s.z),
                if not set(outfilled.getSurrounds(s.x, s.y, s.z)) & set(surrounds_in):
                    # boxels surfacing out is conditionally deleted
                    surface_outs = [e for e in outfilled.getSurrounds(s.x, s.y, s.z) if e.outflag == True]
                    #print "surfacing out: out num = ", len(surface_outs)
                    #for sout in surface_outs:
                    #    print "\t\tout (%d, %d, %d)" % (sout.x, sout.y, sout.z)
                    if len(surface_outs) < 3:
                        newmodel[(s.x, s.y, s.z)] = s
                else:
                    #print ""
                    # boxels surfacing inside always remain
                    newmodel[(s.x, s.y, s.z)] = s
        return newmodel

    #==== Scaled model generation ====#
    def makeMovedModel(self, mvx, mvy, mvz):
        newmodel = BED3D()
        for (x,y,z), b in self.items():
            mb = b.genModPos(x+mvx, y+mvy, z+mvz)
            newmodel[mb.getPosTuple()] = mb
        return newmodel

    #==== rotated model generation ====#
    def makeRotatedModel(self, rotaxis, rotdegree):
        newmodel = BED3D()
        for (x,y,z), b in self.items():
            if rotaxis == 0: # x
                if rotdegree == 0: # 90
                    newx = x
                    newy = -z
                    newz =  y
                elif rotdegree == 1: # -90
                    newx = x
                    newy =  z
                    newz = -y
                elif rotdegree == 2: # 180
                    newx = x
                    newy = -y
                    newz = -z
            elif rotaxis == 1: # y
                if rotdegree == 0: # 90
                    newy = y
                    newx = -z
                    newz =  x
                elif rotdegree == 1: # -90
                    newy = y
                    newx =  z
                    newz = -x
                elif rotdegree == 2: # 180
                    newy = y
                    #newx = -y # bug but interesting!!
                    #newz = -x
                    newx = -x
                    newz = -z
            elif rotaxis == 2: # z
                if rotdegree == 0: # 90
                    newz = z
                    newx = -y
                    newy =  x
                elif rotdegree == 1: # -90
                    newz = z
                    newx =  y
                    newy = -x
                elif rotdegree == 2: # 180
                    newz = z
                    newx = -x
                    newy = -y
            mb = b.genModPos(newx, newy, newz)
            newmodel[mb.getPosTuple()] = mb
        return newmodel

import re
class Boxel(object):
    def __init__(self, x,y,z, r,g,b, alpha):
        self.x = int(x)
        self.y = int(y)
        self.z = int(z)
        #self.color = color
        self.r = int(r)
        self.g = int(g)
        self.b = int(b)
        self.alpha = int(alpha)
        #self.delflag = False
        self.outflag = None # None means not set
        self.doneflag = False

    def getline(self):
        return "%d %d %d %d %d %d %d" % (self.x,self.y,self.z, self.r, self.g, self.b, self.alpha)
    def getPosTuple(self):
        return (self.x, self.y, self.z)

    def __repr__(self):
        return self.getline()

    def genModPos(self,x,y,z):
        newbed = Boxel(x, y, z, self.r, self.g, self.b, self.alpha)
        return newbed

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
    parser = OptionParser(usage=None, version=version) # usageの %prog はOptionParserによってos.path.basename(sys.argv[0])に置換えられる
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
                      default=0,
                      help="specify processing mode: 0(hollow), 1(fill), 2(mirror), 3(scale), 4(move), 5(rotate)",
                      metavar="MODE")
    parser.add_option("-l", "--level",
                      dest="level",
                      type="int",
                      default=1,
                      help="specify level for hollowing mode: should be greater than or equal to 0",
                      metavar="LEVEL")
    parser.add_option("-r", "--mirrormode",
                      dest="mirrormode",
                      type="int",
                      default=0,
                      help="specify mirrormode: 0(use left half), 1(use right half), 2(make other one)",
                      metavar="LEVEL")
    parser.add_option("-s", "--scale",
                      dest="scale",
                      type="int",
                      default=1,
                      help="specify scaling factor integer",
                      metavar="SCALE")
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
    parser.add_option("-a", "--rotaxis",
                      dest="rotaxis",
                      type="int",
                      default=0,
                      help="specify rotation axis: 0(x), 1(y), 2(z)",
                      metavar="ROTAXIS")
    parser.add_option("-d", "--rotdegree",
                      dest="rotdegree",
                      type="int",
                      default=0,
                      help="specify rotation degree: 0(90), 1(-90), 2(180)",
                      metavar="ROTDEGREE")
    (options, args) = parser.parse_args() #引数パーズ

    if options.inputfile == None: 
        parser.print_help()
        exit()
    else:
        return options


if __name__ == "__main__":
    #filename, outname, mode, level, mirrormode, scale = getopt()
    options = getopt()

    bed3D = bed_read(options.inputfile)

    outstr = None
    mode = options.mode
    if mode == 0: # hollow mode
        pruned = bed3D.delClosed(options.level)
        outstr = pruned.printAll()
    elif mode == 1: # fill mode
        filled = bed3D.fillClosed()
        outstr = filled.printAll()
    elif mode == 2: # mirror mode
        newmodel = bed3D.makeMirrorModel(options.mirrormode)
        outstr = newmodel.printAll()
    elif mode == 3: # scaling mode
        newmodel = bed3D.makeScaledModel(options.scale)
        outstr = newmodel.printAll()
    elif mode == 4: # move mode
        newmodel = bed3D.makeMovedModel(options.mvx, options.mvy, options.mvz)
        outstr = newmodel.printAll()
    elif mode == 5: # rotate mode
        newmodel = bed3D.makeRotatedModel(options.rotaxis, options.rotdegree)
        outstr = newmodel.printAll()

    if options.outputfile: 
        with open(options.outputfile, "wt") as of:
            of.write(outstr)
    else:
        import sys
        sys.stdout.write(outstr)

