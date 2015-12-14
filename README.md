
Reproducing the numerical results of the paper
==============================================

The algorithm is described in the academic paper 
[Rigorous enclosures of minimal detectable differences for general ANOVA models](http://reliablecomputing.eu/baharev_anova_power_f_test.pdf).
To reproduce the numerical results of this paper on Linux run:

    mindiffver  <input.txt  >output.txt
    
or on Windows:

    mindiffver.exe  <input.txt  >output.txt

The corresponding executables and the input.txt file are in the 
binary_distribution directory. The executables should work on *any* i386 
compatible architecture. The downside is that they are roughly 7x slower than 
an executable optimized specifically for your CPU.

**The only relevant source file is the `main.cpp` file.** The other 
source files belong to the third party interval arithmetic library C-XSC (see 
the paper).


Example
=======

This program reads from the standard input, and writes to the standard output,
line-by-line. If you have a table of the data to verify in a text file called 
input.txt, then you can run this program like this:

    mindiffver  <input.txt  >output.txt
    
or on Windows:

    mindiffver.exe  <input.txt  >output.txt

The verified results will be written to the output.txt file.

The format of the input.txt and output.txt are detailed right below. 


Input format
============

Each line of the input is supposed to have the following format:

    a  b  x_0  lambda_0  alpha  beta  eps_x  eps_lambda

where the items are separated by arbitrary whitespace, a and b are the shape 
parameters of the noncentral beta distribution, x_0 is the upper alpha quantile
of the (central) beta distribution, beta is the Type II error beta (not 
detecting an effect; power=1-beta), eps_x and eps_lambda are the inflation 
parameters. The search intervals for the correct value of x and lambda are:

    x = [(1-eps_x)*x_0, (1+eps_x)*x_0], and
    lambda = [(1-eps_lambda)*lambda_0, (1+eps_lambda)*lambda_0]. 

Ranges (checked by the program)
------------------------------- 

 -  All input values must be strictly positive
 -  b must be integer
 -  x, alpha, beta, eps_x, eps_lambda < 1 must hold
 -  eps_x >= tol_x, and eps_lambda >= tol_lambda must hold, where 
    tol_x = 10^-12 and tol_lambda = 10^-10 are the currently set tolerances
    for x and lambda in the interval Newton iteration.

Assumptions (not checked)
-------------------------

The parameters are assumed to lie in the domain that is relevant for practical 
applications, roughly: a <= 25, b <= 500, 0.01 <= alpha, beta <= 0.99; the 
inflation parameters are also assumed to be sane, say < 10^-4. Violating these 
assumptions may cause performance degradation and the algorithm may start 
reporting failures but incorrect results will never be produced.


Output format
=============

There are 3 possible outcomes:

a) If the input line contains a solution and the interval Newton method is 
   successful in verifying it, then the output is a line matching the format
   of the input line (items are guaranteed to be tab separated) where x and 
   lambda are guaranteed to have the precision given in the last two 
   columns (currently set to 10^-12 and 10^-10 relative error, respectively).
   
b) If the input line does NOT contain a solution and the interval Newton method
   is successful in verifying it, then the output is a single line saying:
   "The search interval [...] is verified NOT to contain a zero".
   
c) In all other cases a single line error message starting with "Failed ..."
   is printed. See the `failures.txt` input file that systematically triggers 
   all known failure modes, except for the first and the last line of that file
   (those two lines must succeed).


Compiling from source
=====================

You may consider using the binary distribution in which case you avoid the 
hassle of compiling the software. In order to install it from source, 
install C-XSC: 

   http://www.xsc.de/ 

with the default settings, except that both dynamic and static libraries are 
built. **Make sure that all the unit tests of C-XSC pass!** The unit tests are
automatically executed as part of the installation procedure. If you have 
difficulties passing the unit tests, try setting the rounding mode to soft.

Then, in the directory where the source files of the software are, issue the 
following command:

    g++ -O3 -I/path/to/cxsc/include -L/path/to/cxsc/lib *.cpp -Wl,--static -lcxsc -Wl,-Bdynamic -o mindiffver

where the paths /path/to/cxsc/include and /path/to/cxsc/lib are set according
to your C-XSC installation path. This command assumes that you have gcc (g++).

