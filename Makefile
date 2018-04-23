##
 # This file is part of an experimental software implementation of the
 # Erickson-Monma-Veinott algorithm for solving the Steiner problem in graphs.
 # The algorithm runs in edge-linear time and the exponential complexity is
 # restricted to the number of terminal vertices.
 #
 # This version of the source code is released for Track-1 (exact with few
 # terminals) of the Parameterized Algorithms and Computational Experiments
 # Challenge - PACE 2018.  An earlier (multi-threaded) version of this software
 # was released as part of my master's thesis work "Scalable Parameterised
 # Algorithms for two Steiner Problems" at Aalto University, Finland (DOI:
 # http://urn.fi/URN:NBN:fi:aalto-201709046839) and the accompanying software is
 # available in GitHub (https://github.com/suhastheju/steiner-edge-linear).
 #
 # The source code is configured for a gcc build. Other builds are possible but
 # it might require manual configuration of the 'Makefile'.
 #
 # The source code is subject to the following license.
 #
 # Copyright (c) 2018 Suhas Thejaswi
 # Copyright (c) 2017 Suhas Thejaswi
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in all
 # copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 #
 ##

MAKE = make
CC = gcc 
CFLAGS = -O3 -Wall -march=native -std=c99 

SOURCE = reader-el.c

EXE = st-exact

all: $(EXE)

st-exact: $(SOURCE)
	$(CC) $(CFLAGS) -o $@ $< -lm
	
clean:  
	rm -f *.o *.a *~ 
	rm -f $(EXE)
