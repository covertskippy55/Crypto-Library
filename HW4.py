__author__ = 'Dilan'
from math import sqrt


def main():
    modulus = 17  #This is our modulus, default = 17, test mod = 5
    rlist = []  #holds our returnable points
    function = lambda (x): (pow(x, 3, 17) + 4) % modulus  #This is the function for our elliptical curve
    #function = lambda (x): (pow(x,3,5)+(2*x)-1)%modulus #just a test function from the book
    #function = lambda (x): (pow(x,3,5)+(4*x)+4)%modulus #just a test function from the book
    for i in range(0, modulus):  #enumerate all points for the modulus
        squre = sqrt(
            function(i))  # take the square root of the value we got from our function this is because its y^2 =
        #print "Non integer y values: ", squre
        if float.is_integer(squre):  #if the number we received is an integer then we add it to the points we will print
            rlist.append([i, int(squre)])  #add the current x value and the positive y value we got from doing above
            negative = -int(squre)  #This is to take the negative square root value
            if squre != negative:  # if we dont have the same value for negative and positive aka when its 0
                while negative < 0:  #we add modulouses until we reach a positive value
                    negative += modulus
                    negative %= modulus
                rlist.append([i, negative])
    rlist.append([float('inf'),
                  float('inf')])  #add the value of (infinity,infinity) to the list, this is our last point possible
    print rlist  #print our points
    #The output of this file is as follows: [[0, 2], [0, 15], [4, 0], [6, 4], [6, 13], [10, 1], [10, 16], [11, 3], [11, 14], [inf, inf]], where inf indicate infinity
    #These are the enumerated points on the curve y^2 = x^3+4 over modulo 17


main()