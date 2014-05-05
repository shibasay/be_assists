#!/usr/bin/env python
#author  : @shibasay 
#date    : 2014/5/5

header=r"""#be! data file
#(int)x (int)y (int)z (int)r (int)g (int)b (int)a
#delimiters are space' ', tab'\t', camma',' and colon':'."""

class BED3D(dict):
    def printAll(self):
        print header
        for item in self.values():
            print item

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
        xlist_x = sorted([e.x for e in xlist])
        xclosed = xlist_x and len(xlist_x) >= level*2+1 and xlist_x[level-1] < x and xlist_x[-level] > x
        return xclosed
    def checkClosedY(self,x,y,z, level):
        ylist = self.getListZXis(z,x)
        ylist_y = sorted([e.y for e in ylist])
        yclosed = ylist_y and len(ylist_y) >= level*2+1 and ylist_y[level-1] < y and ylist_y[-level] > y
        return yclosed
    def checkClosedZ(self,x,y,z, level):
        zlist = self.getListXYis(x,y)
        zlist_z = sorted([e.z for e in zlist])
        zclosed = zlist_z and len(zlist_z) >= level*2+1 and zlist_z[level-1] < z and zlist_z[-level] > z
        return zclosed
    def checkClosedPos(self,x,y,z, level):
        allclosed = self.checkClosedX(x,y,z, level) and self.checkClosedY(x,y,z, level) and self.checkClosedZ(x,y,z, level)
        return allclosed

    def delClosed(self, level):
        newdata = BED3D()
        for (x,y,z), v in self.items():
            if not self.checkClosedPos(x,y,z, level):
                newdata[(x,y,z)] = v
        return newdata

class BED(object):
    def __init__(self, x, y, z, color, alpha, line):
        self.x = int(x)
        self.y = int(y)
        self.z = int(z)
        self.color = color
        self.alpha = int(alpha)
        self.line = line # original data

    def __repr__(self):
        return self.line

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
                pos, color, alpha = line.split(" ")
                x,y,z = pos.split(",")
                b = BED(x,y,z, color, alpha, line.strip())
                postuple = (int(x), int(y), int(z))
                bed3D[postuple] = b
    return bed3D

if __name__ == "__main__":
    import sys
    filename = sys.argv[1]
    contactlevel = 1
    if len(sys.argv) > 2: 
        contactlevel = int(sys.argv[2])

    bed3D = bed_read(filename)

    pruned = bed3D.delClosed(contactlevel)
    pruned.printAll()

