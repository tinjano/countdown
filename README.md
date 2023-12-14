# countdown
Efficient and unique solutions to Countdown's Numbers game

## Intro
The numbers game on the British TV show [Countdown](https://en.wikipedia.org/wiki/Countdown_(game_show)
involves attempting to combine positive integers into a target number
using addition, positive subtraction, multiplication and division
without remainder. The problem is interesting computationally, even though
it is restricted to exponential complexities, being more complicated than
the subset-sum problem. A naive solution is easy to implement, but is
very inefficient and takes time even for small input sizes. Improvements
can be made by restricting the set of arithmetic expressions being considered,
and some such improvements were suggested in literature (Alliot 2015,
[arXiv](https://arxiv.org/abs/1502.05450), but there is more that can be done.
The rapidly inflating number of possible expressions is the main source
of complexity, so we should aim to exclude all redundancies.

## How it Works
We iterate creating new combinations of expressions but:
- the first operand should be larger than the first, or have a larger address
- for repeated operations, the operands should be decreasing (10+5+2, not 5+10+2)
- no addition after subtraction or multiplication after division
- do not add or subtract an addition or subtraction
- do not multiply or divide with a multiplication or division

There are other possibilities that are technically distinct but are still 
uninteresting and will inflate the set we are working with. Examples are
multiplying/dividing with 1, additions like `a/d+b/d` or `b-a=a`.
Excluding those cases is not enough, as it would still not exclude e.g. `a/d+c+b/d`.
Therefore we can use a dictionary that map masks (vectors that enumerate
how many copies of each input were used; this can be a simple bitmask if
there are not duplicates) to each value that was found. This is the
most expensive part of the program and we may explore the possiblity
of only excluding a few common possibilities like multiplying or dividing by
1, as with smaller input sizes (which are the only option anyway) the
overhead may not be worth it.

## Implementation
A search tree (`struct Map`) is used to keep track of solutions, and 
expressions are represented as binary trees (`struct ExprTree`) so that 
they can be printed in the usual infix form. When iterating through pairs
of expressions, we take care to only produce expressions using the same 
total number of inputs (incl. duplicates). This lets us limit the number
of iterations, and is also conducive to parallelization. SIMD operations
(not implemented) may also be used, but are likely not worth it with 
small sizes.

## How to Use
Firstly, it is not recommended to use input size exceeding 8. If
the memory limit is exceeded (see `MAX_MEMORY_`), the program will exit gracefully.

Compile with 
```
gcc countdown.c -o countdown
```
We input an array of integers &ndash; the first is `N`, the input size
(e.g. 6 numbers in Countdown) followed by `Q`, the number of queries.
The next `N` numbers are the inputs which we can use and the `Q` after that
are the targets we wish to see solution for. Example input and output:
```
./countdown 6 5 100 75 50 25 6 3 952 843 112 19 10998
------
Query: 952
((100+6)*75*3-50)/25
(100+3)*75*6/50+25
2 solutions found.
------
------
Query: 843
(100+50+3)*6-75
(100+25+3)*6+75
2 solutions found.
------
------
Query: 112
100+(75-3)/6
100+50*6/25
100*(25+3)/(75-50)
100*(75+3-50)/25
50*(25+3)*6/75
5 solutions found.
------
------
Query: 19
25-6
100-75-6
75-50-6
75/3-6
25-100*3/50
5 solutions found.
------
------
Query: 10998
0 solutions found.
------
32081 solutions found in total.
Time:0.028443 seconds (28.443000 miliseconds)
Memory freed successfully.
```

## Additional notes
- The speed, depending on the number of solutions, can vary greatly
depending on the inputs. Large, varied and diverse inputs generate more,
while smaller inputs with many duplicates will be much faster with
larger sizes possible. 
- An additional possiblity for optimization is excluding results larger
than a certain threshold. This will exlucde some solutions, but they
may not be interesting to us. With the usual Countdown rules
the largest viable intermediate result is short of 100000.

## Todo
- Comments should be added to the code.

## Known Bugs
- A segmentation fault is possible at the end of the program when freeing
memory. 
