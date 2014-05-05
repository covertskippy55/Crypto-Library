import math
from sys import exit
from random import randrange


# This function will return a list of prime factors using brute force.
def primefactors(n):
    # holds the prime factors
    prime_factors = []
    # we start with the lowest possible factor which is 2
    z = 2
    # we loop until our divisor is bigger than the number we're dividing
    while z <= n:
        # we loop as long as n can be divided without a remainder
        while (n % z) == 0:
            # when we find a number which can divide n without a remainder we add it to our list
            prime_factors.append(z)
            # make n be n divided by z
            n /= z
        # if we dont have a z which divides n perfectly we increment z and try again
        z += 1
    return prime_factors


# This is our Russian peasant based exponentiation function, it takes in the exponent,base and the mod
def exponentiation(x, y, m):
    r = 1  # r starts at 1
    while x != 1:
        # print x,y,r
        # everything here is done mod m allowing for fairly small numbers depending on the size of m
        # if x is divisible by 2 we divide x by 2 and square y
        if x % 2 == 0:
            x = (x / 2)
            y = (y * y) % m
        # otherwise we still halve x and square y but we also mulitply the r and y from the previous row
        else:
            x = (x / 2)
            r = (r * y) % m
            y = (y * y) % m
    ret = (y * r) % m
    # print ("x^y %m =",ret)
    # print "\n"
    # we return when x is 1 and we return the multiplication of y and r mod m
    return ret


# This is our EGCD function
def egcd(x, y):
    q = 0
    u1 = 1  # we start the first row with these values for u1,v1,u2,v2
    v1 = 0
    u2 = 0
    v2 = 1
    first = 1  # This is a flag which allows us not change the Us and the Vs for the first row
    # we go until our remainder is 0
    while y > 0:
        # our check to see if this is the first row if it is dont change Us and Vs
        if first == 1:
            first = 0
        #otherwise we save 2 and v2 in temp variables, then find the u2 and v2 according to the formula,
        # then we set u1 = old u2 and v1 = old v2
        else:
            t1 = u2
            t2 = v2
            u2 = u1 - (t1 * q)
            v2 = v1 - (t2 * q)
            u1 = t1
            v1 = t2
        #the divmod function returns both the quotient and remainder which we store in q and r
        q, r = divmod(x, y)
        #printing of the entire row
        #print ('x = ',x)
        #print ('y = ',y)
        #print ('r = ',r)
        #print ('q = ',q)
        #print ('u1 = ',u1)
        #print ('v1 = ',v1)
        #print ('u2 = ',u2)
        #print ('v2 = ',v2)
        #print "\n"
        #our new x is the old y and the new y is the old r
        x = y
        y = r
    #this is our final answer, which gives us our gcd and u2,v2
    #print ("GCD = ",x)
    #print ("X inverse = ",u2)
    #print ("Y inverse = ",v2)
    #print "\n"
    return x, u2, v2


# This is our factorization by fermats
def fermat_fac(n):
    a = int(math.ceil(math.sqrt(
        n)))  # our a is the ceiling of the sqrt, note that all our calculations are done
    # with whole numbers hence the int conversion.
    b2 = a * a - n  # our b2 is the square of a - the original value
    b = int(
        math.sqrt(b2))  # our b is the square root of b2, this is part of our factor
    # if our b2 is infact a perfect square

    while b * b != b2:  # if b2 is not a perfect square we keep going(if the sqrt of b2 squared is not equal to b2)
        a += 1  # we increment a by 1 and try again for a perfect square
        b2 = a * a - n
        b = int(math.sqrt(b2))
    p = a + b  # out p and q values are given by adding the sqrt of our perfect square to the a value
    # that gave us that perfect square
    q = a - b
    assert n == p * q  # this is basically a test case to make sure out p and q gives us the composite n
    return p, q


# This is the miller rabin method for checking for primes, its slightly different
# than whats taught in class because we only look at one number at a time rather than
# storing them all into an array then looking at the array
def miller_rabin(n, t):
    # The number we're checking the primality has to be bigger than 2
    if n <= 2:
        print "The number you entered is less than 2.\n"
        return
    f = n - 1
    k = 0
    # here we divide f by two repeatedly to find what power of 2 our k is
    while f % 2 == 0:
        k += 1
        f //= 2
    m = f  # our end is whatever is left over from the above loop
    # print ("k = {0}, m = {1}".format(k,m))
    # this is our loop to run it multiple times according to user input, the more times
    # we run the higher the probability
    for _ in range(t):
        a = randrange(1, n - 1)  # generate a random number from 1 to n-1 exclusive
        d, _, _ = egcd(a, n)  # check if our gcd is 1
        if d > 1:
            print "The number is a composite, b/c of gcd"
            return False
        else:
            # if our gcd is 1 we take a^m %n for to get the first element
            x = exponentiation(m, a, n)
            # if the first number is 1 or -1(n-1) we continue
            if x == 1 or x == n - 1:
                continue
            # if the first element was passable we loop k-1 times(because we already did it once)
            for __ in range(k - 1):
                # keep taking the square of our x
                x = exponentiation(2, x, n)
                if x == 1:  #
                    print "This number is a composite b/c during the squaring a non 1 number was found"
                    return False
                # if we find a -1 we can break out of the loop because we already had a 1, so now we have 1 and n-1
                if x == n - 1:
                    break
    print "This number is probably prime"
    return True


# This is our factorization by row method
def rho_factorization(n, s):
    x = y = s  # we start both lists at the given seed
    d = 1
    f = lambda b, m: (b ** 2 + 2) % m  # This is our function which is defined as a lambda
    # we only run until our gcd returns 1 in which case we're done and found a factor
    while d == 1:
        x = f(x, n)  # we do our function to x
        y = f(f(y, n), n)  # and then we do our function twice to y
        d, _, _ = egcd(n, abs(x - y))  # we check the gcd of the difference of x and y with n
        # if the gcd is the same as our number we failed to find a factor of n this time
        if d == n:
            print "Failed to factor, try again with different seed"
            return False
    z = n // d  # this gives us the second factor since rho only finds one by default
    print("The first factor of n is: ", d)
    print ("The second factor of n is: ", z)


# This function will calculate discrete log using baby step giant step method.
# It requires the modulus, a generator of the modulus and the beta.
def baby_giant(p, g, beta):
    m = int(math.ceil(math.sqrt(p - 1)))  # our m is the ceiling of the sqrt of p-1
    print "\n"
    # print m
    # print "\n"
    babystep = []  # holds babystep numbers
    giantstep = []  # holds giantstep numbers
    babystep.append([1,
                     0])  # the first number for babystep is guranteed to be 1 and the index is 0.
    # We store everything as tuples because we need indexes
    f = pow(g, m, p)  # we start by raising g to m mod p
    inv = inverse(f, p)  # this gives us the inverse of the above exponentiation
    #print inverse
    while inv <= 0:  # if the inverse was negative add mods till its positive
        inv += p
    #print inverse
    i = 0
    #we're generating an m number of elements for babystep
    while i <= m:
        #print i
        #print "\n"
        if i != 0:  # this weeds out the first case, because we already did it individually near the top
            babystep.append(
                [pow(g, i, p), i])  # for the rest of it raise g to each i and store that value and the index
        i += 1
    #print babystep
    print "\n"
    #This is our first giant step elemnt, which is just beta at index 0
    giantstep.append([beta, 0])
    i = 0
    # this while loop runs until we reach the end of babystep, this loop is solely to check the
    #  first giant step elemnt against babystep list
    while i < (len(babystep)):
        x = (babystep[i])[0]  # x stores the actual value in the babystep tuple for
        location = (babystep[i])[1]  # location stores the index value in the babystep tuple
        #print x
        #before generating any more elements we check if the first giant step, aka beta, is in the babystep,
        #  if it is we return simply the location
        if x == (giantstep[0])[0]:
            #print (m*0)+location
            return (m * 0) + location
        i += 1
    # now that we checked the first element of giantstep, we generate the rest of the elements and
    #  look through babystep for each element
    i = 0
    while i < m:
        if i != 0:
            y = pow(inverse, i, p)  # exponentiation of the inverse of g^m to i
            y = (beta * y) % p  # multiply beta and y mod p
            giantstep.append([y, i])  # append it first so we can print it even when we find matching elements
            j = 0
            # look through babystep, same loop as above for checking for the first element and
            # if we find it return m*current index + the location from the babystep
            while j < (len(babystep)):
                x = (babystep[j])[0]
                location = (babystep[j])[1]
                #print ("x = ",x)
                #print ("j = ",j)
                #print ("y = ",y)
                if x == y:
                    #print (m*i)+location
                    return (m * i) + location
                j += 1
        i += 1
        #return 0
        #print giantstep


#This is just a short function that uses the egcd function and retrieves just the inverse.
def inverse(x, p):
    ret = egcd(x, p)[1]  # egcd returns multiple things in this order : gcd, inverse of x, inverse of y.
    #  We only want the inverse of x
    if ret < 0:
        ret += p  # if the inverse is negative then add mods until we reach a positive
    return ret


class EllipticalCurve:
    def __init__(self, a, b, p):
        """
        This is our Elliptical Curve class. It will instantiate an elliptical curve of the form x^3+ax+b. It will also
        take a prime if you need to, and it will do all calculations mod said prime. This also has a function parameter
        which is a lambda that will map a given x to the our curve
        :param a: This is our a in ax.
        :param b: This is our b in + b
        :param p: This is our prime number to take mod.
        """
        self.a_input = a
        self.b_input = b
        self.prime = p
        self.function = lambda (x): pow(x, 3, self.prime) + (self.a_input * x) + self.b_input

    def print_curve(self):
        """
        This is a simple print function that takes nothing as input. It will print the curve as a string.

        """
        print "x^3+{0:d}x+{1:d}".format(self.a_input, self.b_input)

    def evaluvate_point(self, x):
        """
        This function will take the given x-coordinate and plug it into the current curve. It will return the y^2 point.
        :param x: This is the x coordinate which we will use to evaluate it on the curve
        :return: The y coordinate soured.
        """
        if self.prime >= 1:  # this if statement is there to ensure that if the prime the user entered is 0,
            # meaning the user wants the curve over the real number field, we wont actually use the modulus

            return self.function(x) % self.prime
        else:
            return self.function(x)

    def point_addition(self, xp, yp, xq, yq):
        """
        This is out point addition function. It uses all the inputs it get in order to calculate the addition of
         two points on the given curve. Not again that if the prime user entered is greater than zero we will take
         everything mod that number.
        :param xp: This is our first x coordinate
        :param yp: This is our first y coordinate
        :param xq: This is our second x coordinate
        :param yq: This is our second y coordinate
        :return: This will return the x and y coordinates of the result of adding the previous two points on an
        the given elliptical curve, or it will return 0 if the modulous was composite and a factor was found
        """
        if self.prime >= 1:
            d, inv, _ = egcd(xq-xp, self.prime)
            if d == 1:
                m = ((yq - yp) * inv) % self.prime  # This will give us the slope of the
                #  line between the two points. Rather than dividing like you would do on a real number field,
                # we must take the inverse of the denominator with respect to the prime and multiply the numerator
                #  with it

                xr = (pow(m, 2, self.prime) - xp - xq) % self.prime
                while xr < 0:
                    xr += self.prime
                yr = m * (xp - xr) - yp % self.prime
                while yr < 0:  # if the y coordinate we got was negative add primes until we get a positive number
                    yr += self.prime
                    yr %= self.prime
            else:
                print "The modulous is not prime, because gcd was > 1."
                r = self.prime/d
                print "One factor of the modulous is: {0:d}, The other factor of the modulous is: {1:d} ".format(d, r)
                return 0
        else:  # This bit is the same as above but without the inverse stuff and the modulus part
            m = ((yq - yp) / (xq - xp))
            xr = pow(m, 2) - xp - xq
            yr = m * (xp - xr) - yp
        return xr, yr

    def evaluate_all_points(self):
        """
        This is our function to evaluate all points on a given elliptical curve. It will start at zero and work until
        we reach our prime - 1. If the curve is not over a prime, there is an infnite amount of possible points and thus
        this cannot be used. In this case it will simply return infnity. This function doesnt check quadratic residues
        so please ignore for now.

        :return: an array of tuples containing the x,y coordinates that belong on the curve
        """
        rlist = []
        for i in range(0, self.prime):
            square_root = math.sqrt(self.evaluvate_point(i))
            if float.is_integer(square_root):  # If the y is not an integer its not on the curve
                rlist.append([i, int(square_root)])
                negative_sqrt = -int(square_root)
                if square_root != negative_sqrt:  # If we have the same y and -y we only need to append it once
                    while negative_sqrt < 0:  # make our negative y into a positive by adding primes
                        negative_sqrt += self.prime
                    negative_sqrt %= self.prime
                    rlist.append([i, negative_sqrt])
        rlist.append([float('inf'), float('inf')])  # The final point in any elliptical curve is infinity, which must be
        # added separately
        return rlist

    def point_doubling(self, xp, yp):
        """
      This function will take the given point on the curve and will add it to itself. We cant do normal point addition
      because the slope will be zero which does not allow for the function to work properly.
      :param xp: The first x coordinate
      :param yp: The first y coordinate
      :return: a tuple that is the new x,y coordinate after point doubling
      """

        if self.prime >= 1:
            d, inv, _ = egcd(2*yp, self.prime)
            if d == 1:
                slope = lambda x, y: ((3 * (pow(x, 2, self.prime)) + self.a_input) * inv) % self.prime
                m = slope(xp, yp)
                xr = (pow(m, 2, self.prime) - 2 * xp) % self.prime
                yr = (m * (xp - xr) - yp) % self.prime
                return xr, yr
            else:
                print "The GCD between the x difference and the modulous is not 1, one factor of modulous is: ", d
                r = self.prime/d
                print "The other factor of the modulous is: ", r
                return 0
        else:
            slope = lambda x, y: (3 * (pow(x, 2)) + self.a_input) / (2 * y)
            m = slope(xp, yp)
            xr = pow(m, 2) - (2 * xp)
            yr = m * (xp - xr) - yp
            return xr, yr


def main():
    """
    This is our main method. It will initialize an infinite loop that will run until the user enters 0 to end the
    program. It will allow the user to chose any of the available cryptographic options. Each loop will also ask the
      user to input any additional information needed to do the function.

    """
    while True:
        choice = int(raw_input(
            "\n Please pick one of the options:\n 0)Exit\n 1)Russian Peasant Algorithm\n 2)Extended GCD\n"
            " 3)Fermat's Factorization\n 4)Miller Rabin Primality Test\n 5)Pollard-Rho Factorization Method\n "
            "6)Baby-Step, Giant-step\n 7)Elliptical Curve Cryptography\n 8)Brute Force Factorization\n:> "))
        if choice == 0:
            exit(0)
        elif choice == 1:
            x = int(raw_input("Please enter the base: "))
            y = int(raw_input("\nPlease enter the exponent: "))
            m = int(raw_input("Please enter the modulus: "))
            ret = exponentiation(y, x, m)
            print ("x^y %m =", ret)
            print "\n"
        elif choice == 2:
            x = int(raw_input("Please enter the x value: "))
            y = int(raw_input("Please enter the y value: "))
            x, u2, v2 = egcd(x, y)
            print ("GCD = ", x)
            print ("X inverse = ", u2)
            print ("Y inverse = ", v2)
            print "\n"
        elif choice == 3:
            x = int(raw_input("please enter the number to be factored: "))
            p, q = fermat_fac(x)
            print("p =", p)
            print("q = ", q)
        elif choice == 4:
            x = int(raw_input("Please enter the number to be checked: "))
            y = int(raw_input(
                "Please enter the number of iterations to do this(The higher this number, "
                "the more accurate your result is): "))
            miller_rabin(x, y)
        elif choice == 5:
            x = int(raw_input("Please enter the number to be checked: "))
            y = int(raw_input("Please enter the seed to be used: "))
            rho_factorization(x, y)
        elif choice == 6:
            x = int(raw_input("Please enter the modulous: "))
            y = int(raw_input("Please enter the generator: "))
            z = int(raw_input("Please enter the beta: "))
            print baby_giant(x, y, z)
        elif choice == 7:
            a = int(raw_input("Please enter the a value: "))
            b = int(raw_input("Please enter the b value: "))
            p = int(raw_input("Please enter the prime your curve is over(if its over Real values, enter 0): "))
            print '\n'
            ecurve = EllipticalCurve(a, b, p)
            while True:
                x = int(raw_input(
                    "What would you like to do? Select the corresponding number\n0) Return to previous menu\n1) "
                    "Print the curve\n2) Evaluate a point on the curve(This will only give you the positive y value, "
                    "the negative is just the negative of it) \n3) Point addition with two points(you MUST have "
                    "two points, if you only have one please select Point Doubling)\n4) Point Doubling(Use this if you "
                    "have only one available point in your curve and you wish to get 2P\n5) Evaluate all points on the "
                    "curve (note that if your prime is very large, this function will fail and return an error,and if "
                    "your curve is not over a prime, this will not work because there are infinite possible points)\n"))
                if x == 0:
                    break
                elif x == 1:
                    ecurve.print_curve()
                elif x == 2:
                    x = int(raw_input("Please enter the x point: "))
                    y2 = ecurve.evaluvate_point(x)  # The return value we get evaluvate_point is actually y^2, so we
                    # take the squire root of it. If its an integer then we have a point on the curve, so output the x
                    # and y coordinates. Otherwise let the user know its not on the curve.
                    if float.is_integer(math.sqrt(y2)):
                        print [x, int(math.sqrt(y2))]
                    else:
                        print 'The x coordinate you entered is not within the curve\n'
                elif x == 3:
                    xp = int(raw_input("Please enter the first x coordinate: "))
                    yp = int(raw_input("Please enter the first y coordinate: "))
                    xq = int(raw_input("Please enter the second x coordinate: "))
                    yq = int(raw_input("Plesae enter the second y coordinate: "))
                    print ecurve.point_addition(xp, yp, xq, yq)
                elif x == 4:
                    xp = int(raw_input("Please enter the x - coordinate: "))
                    yp = int(raw_input("Please enter the y - coordinate: "))
                    print ecurve.point_doubling(xp, yp)
                elif x == 5:
                   # print ecurve.evaluate_all_points()
                    print "This function is currently being reconstructed."
        elif choice == 8:
            x = int(raw_input("Please enter the number to be factored: "))
            print primefactors(x)


if __name__ == '__main__':  # this is the boilerplate portion
    main()