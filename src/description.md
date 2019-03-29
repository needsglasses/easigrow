# easiGRO

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

```
easigro -q ft55.seq -s 300 -r -b seft-newman84 --cycle_method rainflow -d walker:default -a 10e-6 -e 5e-3 -o block,a -n 10 
```

This means:

Read in the sequence with the filename 'ft55.seq' [-q ft55.seq]
scaling by a factor of 300 MPa [-s 300] and reorder the sequence [-r]
to close all cycles. Use the beta model for a 'semi-elliptical surface
crack in a finite plate in tension' by Newman and Raju [-b
seft-newman84] and calculate the crack size by summing up the growth
increment for each rainflow cycle [--cycle_method rainflow] using the
Walker da/dN equation with default material parameters [-d
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

