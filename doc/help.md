This is the output from `easigrow --help`.

```
easiGrow: A crack growth modelling program 0.0.0
Paul White

# easigrow

     Calculates fatigue crack growth rates and finds the
     optimum parameters of a crack growth model to match predictions
     with measurements.

     The distances to the free edges on the coupon are specified
     relative to the origin of the crack. i.e. they are not the coupon
     dimensions. The crack depth 'a' is in the `forward`
     direction. The crack length 'c' is in the `sideways`
     direction. If the crack is not symmetrically located within the
     coupon, use the distance to the nearest free edge.

     a)        |<-2c->|               b)       |<c>|                         
            ________________  ___               ________________  ___     
           |    \     / ^   |  ^               |   |   ^        |  ^      
           |     \___/  v a |  forward         |_ /    v a      |  forward
           |                |                  |                |         
           |________________| _v_              |________________| _v_     

                   |sideways|                  |<-- sideways -->|         

     Figure: Cross section of a coupon showing crtical dimension names
     with a) semi-eliptical crack, and b) corner crack.
	 
Example
========

easigrow -q ft55.seq -s 300 -r -b seft-newman84 --cycle_method rainflow -d walker:default -a 10e-6 -e 5e-3 -o block,a -n
10 


This means:

Read in the sequence with the filename 'ft55.seq' [-q ft55.seq]
scaling by a factor of 300 MPa [-s 300] and reorder the sequence [-r]
to close all cycles. Use the beta model for a 'semi-elliptical surface
crack in a finite plate in tension' by Newman and Raju [-b
seft-newman84] and calculate the crack size by summing up the growth
increment for each rainflow cycle [--cycle_method rainflow] using the
Walker da/dN equation with default material parameters [-m
walker:default].  Starting at an initial crack size 10 um [-a 10e-6]
grow the crack until a final size 5 mm [-e 5e-3], writing out the
variables 'block' and 'a' [-o block,a] every 10 blocks [-n 10].

# How the program works

Think of the program flow as

1. Read in data
2. Filter the sequence (turning point, rainflow, risefall, deadband etc.) and convert to cycles
3. Filter the list of cycles
4. If required, optimise any matching parameters
5. Perform a crack growth calculation
6. Write out requested output


USAGE:
    easigrow [FLAGS] [OPTIONS]

FLAGS:

    -h, --help             Prints help information
    -l, --list             list all available output variables, growth, beta and da/dN models
        --seq_cyclemods    only keep those turning points used in remaining cycles
    -r, --seq_reorder      re-order the sequence to close cycles (default: no reordering)
        --seq_tp           remove non-turning points from the sequence (default: true)
        --summary          print summary of sequence
    -V, --version          Prints version information
    -v, --verbose          print more information

OPTIONS:
    -e, --limit_a <LENGTH>             final termination crack sizes (default 1e-3 m)
    -a, --astart <LENGTH>              initial crack sizes (default 10e-6 m)
    -b, --beta <NAME>                  beta geometry model (default seft-newman84)
        --beta_outfile <FILE>          write the beta model to a file
    -c, --crack_infile <FILE>          crack growth measurements that will be matched during optimisation
        --crack_weight <WEIGHT>        weighting factor for matching crack growth curve (default 1.0)
        --cycle_deadband <MIN,MAX>     eliminates cycles with turning points within this band
        --cycle_infile <FILE>          read the cycles from an AFGROW input file
        --cycle_max <MAX>              any range greater than this will be set to this value
        --cycle_method <METHOD>        method used to extract cycles from sequence (default rainflow) [possible values:
                                       rainflow, tension]
        --cycle_min <MIN>              any range less than this will be set to this
        --cycle_outfile <FILE>         write the modified sequence to a file
        --cycle_rem_big <DELTA>        remove cycles with a range bigger than this
        --cycle_rem_small <DELTA>      remove cycles with a range smaller than this
    -d, --dadn <NAME>                  dadn material equation (default white:barter14_aa7050-t7451)
        --forward <DISTANCE>           distance forwards from crack origin to free edge in 'a' direction of growth
                                       (default: infinite)
        --image_bar <LENGTH>           size of scale bar for image (default 50e-6 m)
        --image_outfile <FILE>         generate a pseudo image of the fracture surface and write to FILE
        --image_size <V,H>             size of image VxH pixels (default 300x8000)
        --image_type <TYPE>            type of image output (default sem) [possible values: Sem, Optical]
    -N, --limit_block <N>              maximum number of blocks (default +1000)
        --limit_k1c <KIC>              fracture toughness K1C for termination (default 33 MPa sqrt(m))
        --limit_yield <STRESS>         yield or flow stress for plastic zone calculations and net section yield (default
                                       450 MPa).
        --opt_infile <FILE>            optimises model parameters to match a set of target crack growth data.
        --opt_max <N>                  maximum number of iterations to be used for optimisation (default 100)
        --opt_method <METHOD>          method used for optimisation (default Levenberg) [possible values: Sweep, Nelder,
                                       All]
        --opt_nelder <p1,...,p5>       nelder-mead parameters (default step: 0.1, alpha: 1.0, gamma: 2.0, rho: 0.5,
                                       sigma: 0.5)
        --opt_sweep <f1,...,fM>        perform a brute force sweep over a range of scaling factors applied to M dadn
                                       parameters
        --opt_tol <TOL>                terminate optimisation when the change in error function is less than this amount
                                       (default 1e-3)
    -n, --output_every <N>             output crack growth data every +N blocks or -N cycles (default +1)
        --output_lines <l1,l2,...>     output growth at the given load lines l1,l2,l3... (default 1)
    -o, --output_vars <v1,v2,...>      comma separated list of variables to be output (default block,a)
    -p, --parameters <p1,p2,...,pM>    material parameters (default parameters for white:barter14_aa7050-t7451)
        --radius <R>                   radius of hole or notch (default: infinite)
    -s, --scale <SCALE>                scale factor for load sequence (no default)
    -q, --seq_infile <FILE>            file of loading sequence for one block of loading
        --seq_max <MAX>                any value greater than this will be set to this value
        --seq_min <MIN>                any value less than this will be set to this
        --seq_outfile <FILE>           write the modified sequence to a file
        --seq_rem_big <BIG>            any values greater than this will be removed
        --seq_rem_small <SMALL>        any value less than this will be removed
        --sequence <S1,S2,...>         explicitly set the sequence of loading (default: [0.0, 1.0, 0.0])
        --sideways <DISTANCE>          distance sideways from crack origin to free edge in 'c' direction of growth
                                       (default: infinite)
        --youngs <MODULUS>             Young's modulus for compact tension coupon displacements (default 71000 MPa)
```
