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
#        top = self.get((x,y,z), NONE)
#        bot = self.get((x,y,z), NONE)
#        lef = self.get((x+1,y,z), NONE)
#        rig = self.get((x-1,y,z), NONE)
#        fore = self.get((x,y,z+1), NONE)
#        back = self.get((x,y,z-1), NONE)
#        return top,bot,lef,rig,fore,back

    def checkClosedX(self,x,y,z,  level):
        xlist = self.getListYZis(y,z)
        xlist_smaller = sorted([e.x for e in xlist if e.x < x])
        xlist_larger  = sorted([e.x for e in xlist if e.x > x])
        xclosed = len(xlist_smaller) >= level and len(xlist_larger) >= level
        return xclosed
    def checkClosedY(self,x,y,z, level):
        ylist = self.getListZXis(z,x)
        ylist_smaller = sorted([e.y for e in ylist if e.y < y])
        ylist_larger  = sorted([e.y for e in ylist if e.y > y])
        yclosed = len(ylist_smaller) >= level and len(ylist_larger) >= level
        return yclosed
    def checkClosedZ(self,x,y,z, level):
        zlist = self.getListXYis(x,y)
        zlist_smaller = sorted([e.z for e in zlist if e.z < z])
        zlist_larger  = sorted([e.z for e in zlist if e.z > z])
        zclosed = len(zlist_smaller) >= level and len(zlist_larger) >= level
        return zclosed
    def checkClosedPos(self,x,y,z, level):
        allclosed = self.checkClosedX(x,y,z, level) and self.checkClosedY(x,y,z, level) and self.checkClosedZ(x,y,z, level)
        return allclosed

    def delClosed(self, level):
        newdata = BED3D()
        for (x,y,z), v in self.items():
            if level == 0 or (not self.checkClosedPos(x,y,z, level)):
                newdata[(x,y,z)] = v
        return newdata

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
                    elif self.checkClosedPos(x,y,z, 1):
                        fillboxel = "%d,%d,%d %d,%d,%d %d" % (x,y,z, FILL_R, FILL_G, FILL_B, FILL_A)
                        newdata[(x,y,z)] = BED(fillboxel)
                        #print "fill!" 
                    else:
                        #print "do nothing..."
                        pass
        return newdata

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
        print "@delHalf: width = %d, offset = %d, center = %d, xmin = %d, xmax = %d" % (width, offset, center, self.xmin, self.xmax)
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
        print "@makeMirror: width = %d, oddflag = %d, xmin = %d, xmax = %d" % (width, oddflag, self.xmin, self.xmax)
        for (x,y,z), v in self.items():
            newdata[(x,y,z)] = v
            diff = x-center
            newx = center-diff+offset if lr==0 else center-diff-offset
            if newx != x:
                newdata[(newx,y,z)] = v.genModPos(newx,y,z)
        return newdata

class BED(object):
    def __init__(self, line):
        pos, color, alpha = line.split(" ")
        x,y,z = pos.split(",")
        self.x = int(x)
        self.y = int(y)
        self.z = int(z)
        self.color = color
        self.alpha = int(alpha)
        self.line = line # original data

    def __repr__(self):
        return self.line

    def genModPos(self,x,y,z):
        newbed = BED(self.line)
        newbed.x = x
        newbed.y = y
        newbed.z = z
        newbed.line = "%d,%d,%d %s %d" % (x,y,z, self.color, self.alpha)
        return newbed

#    def is_top_of(self, other):
#        return self.y == other.y-1
#    def is_bottom_of(self, other):
#        return self.y == other.y+1
#
#    def is_left_of(self, other):
#        return self.x == other.x-1
#    def is_right_of(self, other):
#        return self.x == other.x+1
#
#    def is_fore_of(self, other): # foreground
#        return self.z == other.z-1
#    def is_back_of(self, other): # background
#        return self.z == other.z+1

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
        pruned = bed3D.delClosed(level)
        outstr = pruned.printAll()
    elif mode == 1: # fill mode
        filled = bed3D.fillClosed()
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

