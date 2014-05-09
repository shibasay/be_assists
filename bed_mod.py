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
            outlist.append(v.line)
        return "\n".join(outlist)

#    def getListXis(self, xpos):
#        return [v for (x,y,z), v in self.items() if x==xpos]
#    def getListYis(self, ypos):
#        return [v for (x,y,z), v in self.items() if y==ypos]
#    def getListZis(self, zpos):
#        return [v for (x,y,z), v in self.items() if z==zpos]

    def getListXYis(self, xpos, ypos):
        return [v for (x,y,z), v in self.items() if x==xpos and y==ypos]
    def getListYZis(self, ypos, zpos):
        return [v for (x,y,z), v in self.items() if y==ypos and z==zpos]
    def getListZXis(self, zpos, xpos):
        return [v for (x,y,z), v in self.items() if z==zpos and x==xpos]

#    def getSurrounds(self, x,y,z): 
#        top = self.get((x,y,z), None)
#        bot = self.get((x,y,z), None)
#        lef = self.get((x+1,y,z), None)
#        rig = self.get((x-1,y,z), None)
#        fore = self.get((x,y,z+1), None)
#        back = self.get((x,y,z-1), None)
#        return top,bot,lef,rig,fore,back

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

    #============= adjucent existence based delete decision ==========#
    def checkClosed_base(self, value, level, check_strong, smaller_list, larger_list):
        def fand(x,y): return (x and y)
        closed_weak = len(smaller_list) >= level and len(larger_list) >= level
        smaller_adj_exist = True
        larger_adj_exist = True
        if check_strong and closed_weak:
            smaller_adj_exist = reduce(fand, [smaller_list[-i] == value-i for i in range(1, level+1)], True)
            larger_adj_exist = reduce(fand, [larger_list[i-1] == value+i for i in range(1, level+1)], True)
        closed = closed_weak and smaller_adj_exist and larger_adj_exist
        return closed
    def checkClosedX(self,x,y,z,  level, check_strong):
        xlist = self.getListYZis(y,z)
        xlist_smaller = sorted([e.x for e in xlist if e.x < x])
        xlist_larger  = sorted([e.x for e in xlist if e.x > x])
        return self.checkClosed_base(x, level, check_strong, xlist_smaller, xlist_larger)
    def checkClosedY(self,x,y,z, level, check_strong):
        ylist = self.getListZXis(z,x)
        ylist_smaller = sorted([e.y for e in ylist if e.y < y])
        ylist_larger  = sorted([e.y for e in ylist if e.y > y])
        return self.checkClosed_base(y, level, check_strong, ylist_smaller, ylist_larger)
    def checkClosedZ(self,x,y,z, level, check_strong):
        zlist = self.getListXYis(x,y)
        zlist_smaller = sorted([e.z for e in zlist if e.z < z])
        zlist_larger  = sorted([e.z for e in zlist if e.z > z])
        return self.checkClosed_base(z, level, check_strong, zlist_smaller, zlist_larger)
    def checkClosedPos(self,x,y,z, level, check_strong):
        allclosed = (self.checkClosedX(x,y,z, level, check_strong) and
                     self.checkClosedY(x,y,z, level, check_strong) and
                     self.checkClosedZ(x,y,z, level, check_strong))
        return allclosed

    def setDelFlag(self, level):
        for (x,y,z), v in self.items():
            if level > 0 and self.checkClosedPos(x,y,z, level, True):
                v.delflag = True

    def delClosed(self, level):
        self.setDelFlag(level)
        newdata = BED3D()
        for (x,y,z), v in self.items():
            if v.delflag == False:
                newdata[(x,y,z)] = v
        return newdata

    def fillClosed(self):
        self.updatePosMaxMin()
        newdata = BED3D()
        for x in range(self.xmin, self.xmax+1):
            for y in range(self.ymin, self.ymax+1):
                for z in range(self.zmin, self.zmax+1):
                    #print "check (%d, %d, %d)" % (x,y,z), 
                    v = self.get((x,y,z), None)
                    if v: 
                        newdata[(x,y,z)] = v
                        #print "exist original data" 
                    elif self.checkClosedPos(x,y,z, 1, False):
                        fillboxel = "%d,%d,%d %d,%d,%d %d" % (x,y,z, FILL_R, FILL_G, FILL_B, FILL_A)
                        newdata[(x,y,z)] = BED(fillboxel)
                        #print "fill!" 
                    else:
                        #print "do nothing..."
                        pass
        return newdata

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
        xlist_smaller = sorted([(e.x, e) for e in xlist if e.x < x])
        xlist_larger  = sorted([(e.x, e) for e in xlist if e.x > x])
        return self.checkClosed_base_usingOut(x, level, xlist_smaller, xlist_larger)
    def checkClosedY_usingOut(self,x,y,z, level):
        ylist = self.getListZXis(z,x)
        ylist_smaller = sorted([(e.y, e) for e in ylist if e.y < y])
        ylist_larger  = sorted([(e.y, e) for e in ylist if e.y > y])
        return self.checkClosed_base_usingOut(y, level, ylist_smaller, ylist_larger)
    def checkClosedZ_usingOut(self,x,y,z, level):
        zlist = self.getListXYis(x,y)
        zlist_smaller = sorted([(e.z, e) for e in zlist if e.z < z])
        zlist_larger  = sorted([(e.z, e) for e in zlist if e.z > z])
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
        for k, v in self.items():
            if v.outflag == True:
                del(self[k])

    def getFillingOut(self,x,y,z):
        if (self.xmin <= x <= self.xmax and 
            self.ymin <= y <= self.ymax and
            self.zmin <= z <= self.zmax):
            b = self.get((x,y,z), None)
            if b == None:
                newline = "%d %d %d %d %d %d %d" % (x, y, z, FILL_R, FILL_G, FILL_B, FILL_A)
                b = BED(newline)
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
    def setOutBEDs(self):
        self.initOutFlag()
        self.updatePosMaxMin()
        # make base position which is absolutely OUT
        newline = "%d %d %d %d %d %d %d" % (self.xmin-1, self.ymin-1, self.zmin-1, 
                                            FILL_R, FILL_G, FILL_B, FILL_A)
        n = BED(newline)
        n.outflag = True
        self[(self.xmin-1, self.ymin-1, self.zmin-1)] = n
        self.updatePosMaxMin()

        undones = self.getSurroundingUndoneOuts(n)
        while undones:
            newundones = set()
            for b in undones:
                for e in self.getSurroundingUndoneOuts(b):
                    newundones.add(e)
            undones = newundones

    def delClosed_usingOut(self, level):
        newdata = BED3D()
        self.setOutBEDs()
        for (x,y,z), v in self.items():
            if v.outflag != True and not self.checkClosedPos_usingOut(x,y,z, level):
            #if True: # test
                newdata[(x,y,z)] = v
        return newdata

    def fillClosed_usingOut(self):
        self.setOutBEDs()
        self.updatePosMaxMin()
        newdata = BED3D()
        for x in range(self.xmin, self.xmax+1):
            for y in range(self.ymin, self.ymax+1):
                for z in range(self.zmin, self.zmax+1):
                    #print "check (%d, %d, %d)" % (x,y,z), 
                    v = self.get((x,y,z), None)
                    if v:
                        if not v.outflag: 
                            newdata[(x,y,z)] = v
                            #print "exist original data" 
                    else:
                        fillboxel = "%d,%d,%d %d,%d,%d %d" % (x,y,z, FILL_R, FILL_G, FILL_B, FILL_A)
                        newdata[(x,y,z)] = BED(fillboxel)
                        #print "fill!" 
        return newdata

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

import re
class BED(object):
    def __init__(self, line):
        x,y,z, r,g,b, alpha = re.split(r"[ ,;\t]", line)
        self.x = int(x)
        self.y = int(y)
        self.z = int(z)
        #self.color = color
        self.r = int(r)
        self.g = int(g)
        self.b = int(b)
        self.alpha = int(alpha)
        self.line = line # original data
        self.delflag = False
        self.outflag = None # None means not set
        self.doneflag = False

    def __repr__(self):
        return self.line

    def genModPos(self,x,y,z):
        newbed = BED(self.line)
        newbed.x = x
        newbed.y = y
        newbed.z = z
        newbed.line = "%d %d %d %d %d %d %d" % (x,y,z, self.r, self.g, self.b, self.alpha)
        return newbed

def bed_read(bedfilename):
    bed3D = BED3D()
    with open(bedfilename) as f:
        for line in f:
            if line[0] != '#':
                b = BED(line.strip())
                postuple = (b.x, b.y, b.z)
                bed3D[postuple] = b
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
                      help="specify processing mode: 0(hollow) or 1(fill)",
                      metavar="MODE")
    parser.add_option("-l", "--level",
                      dest="level",
                      help="specify level for hollowing mode: should be greater than or equal to 0",
                      metavar="LEVEL")
    parser.add_option("-r", "--mirrormode",
                      dest="mirrormode",
                      help="specify mirrormode 0: use left half, 1: use right half, 2: make other one",
                      metavar="LEVEL")
    (options, args) = parser.parse_args() #引数パーズ

    infile  = options.inputfile  if options.inputfile else None
    outfile = options.outputfile if options.outputfile else None
    mode    = int(options.mode)  if options.mode else 0
    level   = int(options.level) if options.level else 1
    mirrormode = int(options.mirrormode) if options.mirrormode else 0

    if infile == None: 
        parser.print_help()
        exit()
    else:
        return (infile, outfile, mode, level, mirrormode)


if __name__ == "__main__":
    filename, outname, mode, level, mirrormode = getopt()

    bed3D = bed_read(filename)

    outstr = None
    if mode == 0: # hollow mode
        pruned = bed3D.delClosed_usingOut(level)
        outstr = pruned.printAll()
    elif mode == 1: # fill mode
        filled = bed3D.fillClosed_usingOut()
        outstr = filled.printAll()
    elif mode== 2: # mirror mode
        newmodel = bed3D.makeMirrorModel(mirrormode)
        outstr = newmodel.printAll()

    if outname: 
        with open(outname, "wt") as of:
            of.write(outstr)
    else:
        import sys
        sys.stdout.write(outstr)

