Overview
--------
This software repository contains an experimental software implementation of the
Erickson-Monma-Veinott algorithm for solving the Steiner problem in graphs. It
is a parameterised algorithm, which runs in edge-linear time and the exponential
complexity is restricted to the number of terminals. The software is written in
C programming language.

This version of the source code is realeased for Track-1 (exact with few
terminals) of the Parameterized Algorithms and Computational Experiments
Challenge - PACE 2018. An earlier (multi-threaded) version of this software was
released as a part of my master's thesis work "Scalable Parameterised Algorithms
for two Steiner Problems" at Aalto University, Finland (DOI:
http://urn.fi/URN:NBN:fi:aalto-201709046839) and the accompanying software is
available in GitHub (https://github.com/suhastheju/steiner-edge-linear).

The source code is configured for a gcc build. Other builds are possible but it
might require manual configuration of the 'Makefile'.

The source code is subject to MIT license, see 'LICENSE.txt' for details.

Building
--------
Use GNU make to build the software. Check 'Makefile' for more details.

Using the software
------------------
Usage: st-exact -s <seed> <input-graph>

Mandatory arguments:
    seed            integer value in the range 0 to 2^31
    input-graph     input graph in DIMACS STP format

Output:
    VALUE           cost of minimum Steiner tree
    edge-list       edge-list of a minimum Steiner tree

Sample output: st-exact -s 1234 instance001.gr

VALUE 503
47 25
25 1
47 53
53 43
43 22
22 28
28 8
8 29
29 7
7 9
29 17
17 24
24 40
