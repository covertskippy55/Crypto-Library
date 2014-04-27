__author__ = 'Dilan'
from math import sqrt


def main():
    modulous = 5
    rlist = []
    function = lambda (x): (pow(x, 3, 17) + 4) % 17
    #function = lambda (x): (pow(x,3,5)+(4*x)+4) % modulous
    for i in range(0, 16):
    #for i in range(0,5):
        squre = sqrt(function(i))
        #print "Non integer y values: ", squre
        if float.is_integer(squre):
            rlist.append([i, int(squre)])
            negative = -int(squre)
            if negative != squre:
                while negative <=0:
                    negative+=modulous
                    negative %= modulous
                rlist.append([i,negative])
    rlist.append([float('inf'),float('inf')])
    print rlist


main()