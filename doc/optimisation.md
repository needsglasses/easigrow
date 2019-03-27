# Optimisation
To calculate the optimium coefficients that produces crack growth that matches
existing crack growth data we need to supply the following information:

1. a crack growth model corresponding to each crack measurement file
  that including beta, stress, start and end cracks size, dadn model
2. loading sequences used for the crack measurements

Each line in the optimisation file is an easigrow command line that will
generate a crack growth curve along with a measurement file of crack
growth data that it tries to match.

```


```

Here is a sample listing of crack growth measurements obtained from a microscope.

```
Item    X (mm)  Y (mm)   Z (mm)     Radius (mm)
1.00     0.0000     0.0000  -0.0000  0.0000
2.00     0.1280     0.1638   0.1991  0.2079
3.00    -0.0128     6.0006   0.0009  6.0007
4.00    -3.7212     6.2297  -0.0620  7.2565
5.00    -3.7212     6.2297  -0.0620  7.2565
6.00    -3.7572     6.2297  -0.0666  7.2750
7.00    -3.7910     6.2280  -0.0670  7.2911
8.00    -3.6810     5.7518  -0.1105  6.8288
9.00    -3.6810     5.7518  -0.1105  6.8288
10.00   -3.7184     5.7518  -0.0871  6.8490
```

Note, these measurements are taken at the same point within a block or
feature on the fracture surface. In this case, because we are only
dealing with blocks, the exact load line the measurements were made at
is not important. Unless otherwise specified, the program will assume
it is the first load line. Contiguous measurements are made until the
blocks can no longer be identified. We indicate that a measurement is
not contiguous with the previous block and is the start of a new run
by pressing the measurement button twice. Here measurements 4 and 5,
and 8 and 9 are both the same indicating the start of a new run.

This measurement file is turned into an input file that looks like
```bash
0.2079

6.0007

7.2565
7.2750
7.2911

6.8288
6.8490
```

An example of the output from an optimisation run is shown below. The
output shows how the optimisation is doing. Here we have used the
Walker equation which contains three variables which are optimised,
they are [0.00000001, 3, 0.5]. We also include the crack growth from 4
coupons, each with an associated crack growth file. The error term is
calculated for each of the files and compared with the
easigrow predictions using the initial set of parameters. In this
case the initial parameters are the default values which are not very
good.

There are a number of different optimisation methds but the simplest
is Nelder optimsiation method which is used for this example. At the
start of a Nelder optimisation the algorithm adds a small increment to
each of the optimisation parameters in turn to determine which direction is
best to begin the optimisation. Then it will head off in
the most favourable direction adjusting the other parameters looking
for directions that result in some reduction in the overall
error. When the level of improvement in each step falls below the
convergence threshold it will terminate the optimisation with a message saying it has
converged.

The program may also terminate if it has exceeded the maximum number
of iterations. The default value for the maximum number of iterations
in an optimisation is 100. This can be overridden on the command line,
as we have done here, increasing the maximum iteration limit to
1000. The initial error terms for each of the crack growth files are
all in the region of 100--200 i.e. not a good prediction. Eventually,
the solution converges and we have an *optimised* answer where the
error between the optimised predictions and the measured crack growth
are of the order of 20--30. A significant improvement. The final
optimised parameters are [1.23759e-9, 2.0023, 1.8756].

Inside the program, the optimisation is performed on a
non-dimensionalised version of the parameters, where the parameters
are normalised based on their starting values. This makes the scaling
of the parameters much closer to each other, which is better for most
optimisation routines. The resultant normalised factors are [0.12376,
0.66743, 3.7511].

```
easigro --optimise optim --opt_max 1000 --opt_method Nelder -m walker:default

```

