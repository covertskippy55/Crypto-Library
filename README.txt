This is a list of cryptographic functions written in python for my CSI 426 class at UAlbany. This interactive program allows you to chose from the following working functions:
Prime Factorization- Three different methods are available.
	1)Brute force prime factorization will provide you  with a list of prime factors, however it is the slowest because it simply checks all numbers to see if they are part of the factorization.
	2)Fermat's factorization will allow you to factor a composite number that is the product of two primes that are fairly close to each other.
	3)Rho Factorization uses a probabilistic approach. It uses a user inputted seed. This method can fail in which case you simply enter a new seed and try again. Eventually you will get two factors assuming your number is a factorable number.


	
Exponentiation using Russian Peasant Algorithm, which uses squaring and halving modulo p to achieve exponentiation


Extended Euclidean Greatest Common Denominator Algorithm which will return the greatest common denominator between two numbers as well as the inverse of both with respect to the other. 

Miller Rabin method of primality testing, which uses a probabilistic approach to determine if a number is a prime number. Miller Rabin is particularly good at dealing with Carmichael Numbers.

Discrete Logarithm done using Baby-step, Giant-step algorithm. Pohlig-Hellman is a possibility for a new logarithm function.

An Elliptical Curve class, that will create an elliptical curve given the user input of a, b and a prime. The prime can be 0 if you wish to create the curve over the real number field. The following functions are available once you create your elliptical curve:
1) Print Curve will print the current curve in the form x^3+ax+b
2) Evaluate point will evluvate a point on the curve given an x-coordinate.
3) Point Addition will add two points within the curve. The user must input two points for this to work. If you only have 1 point, you must select point doubling option. 
3) Evaluate all points on the curve will print all of the x,y coordinates within the curve. If your curve is over the real numbers, there is an infnite number of possible points, thus the program will simply print the identity of the curve, infinity. 	
