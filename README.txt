This is a list of cryptographic functions written in python for my CSi 426 class at UAlbany. The interactive program allows you to chose from the following working functions:
Prime Factorization- Three different methods are available.
	1)Brute force prime factorization will provide you  with a list of prime factors, however it is the slowest.
	2)Fermat's factorization will allow you to factor a composite number that is the product of two primes that are fairly close to each other.
	3)Rho Factorization uses a probabilistic uses a user inputted seed. This method can fail in which case you simply enter a new seed and try again. Eventually you will get two factors assuming your number is a factorable number.


	
Exponentiation using Russian Peasant Algorithm, which uses squaring and halving modulo p to achieve exponentiation


Extended Euclidean Greatest Common Denominator Algorithm which will return the greatest common denominator between two numbers as well as the inverse of both with respect to the other. 

Miller Rabin method of primality testing, which uses a probabilistic approach to determine if a number is a prime number. Miller Rabin is particularly good at dealing with Carmichael Numbers.

Discrete Logarithm done using Baby-step, Giant-step algorithm. I am currently working on implementing Pholig-Hellman as an alternate method. Baby-step, Giant-step tends to be space inefficient but given current available space its not a big issue.   	
	
