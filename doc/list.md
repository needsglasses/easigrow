The is the output from `easigrow --list`, slighlty reformatted for markdown.

# Program Highlights

       * **Sequence filtering** Performs sequence reordering,
          turning-point, dead-band, rise-fall filtering and rain-flow
          counting.  The filtered sequences may be written to a file.

       * **Inbuilt data** Comes with a selection of beta factors,
          material data and crack growth models.

       * **Calculated parameters** Calculates additional parameters for
          characterising the state of the crack tip so that better
          crack growth equations can be developed based on the most
          applicable parameters for a material.

       * **Optimisation** Optimises the crack growth model parameters to
          minimise the difference between predicted and measured crack
          growth. The measured data need not be for the entire history
          i.e., one or more fractographic measurements of the width of
          a block.

       * **Image generation** Generates a pseudo fractographic image of the fracture
          surface to see how easy it is to identify the individual
          blocks. 

# Units

  The internal da/dN data are all in units for stress intensity of
  (MPa m^0.5) and growth in (m). Most beta equations use an applied
  far field stress which is in units of (MPa). However, the
  compact-tension beta factor compact_tada73 uses applied load not stress and
  is in units of load of (MN) and coupon dimensions are in (m). The
  width and depth will need to be set for the compact-tension beta
  function otherwise **easiGro** will assume an infinite plate and the
  crack will not grow.

# Output Parameters

block                block number
line                 line number
a/c                  crack aspect ratio
a/d                  cracked fraction of forward distance
c/b                  cracked fraction of sideways distance
k                    current cycle K
dk                   current cycle dK
r                    current cycle R
beta_a               beta factor at a
beta_c               beta factor at c
a                    forward distance of crack from origin
c                    sideways distance of crack from origin
da                   crack growth increment at a
dc                   crack growth increment at c
mono                 largest monotonic plastic zone 
cyclic               largest cyclic plastic zone 
a/mono               ratio of crack length to monotonic plastic zone size
a/cyclic             ratio of crack length to cyclic plastic zone size
mono/da              ratio of current cyclic plastic zone to current da
cyclic/da            ratio of current cyclic plastic zone to current da
peak                 scaled peak stress of current cycle
valley               scaled valley stress of current cycle

# Beta Models

qct-broek86                                         quarter circular crack in an infinite plate in tension [broek86]
seft-newman84        a/d, a/c, c/b, phi             semi-elliptical surface crack in a finite plate in tension [Newman79]
seit-anderson05      a/c, phi                       semi-elliptical surface crack in an infinite plate in tension [Anderson05]
qcft-murakami87      a/d                            quarter circular corner crack in a finite plate in tension [Murakami87]
qeft-newman84        a/d, a/c, c/b, phi             quarter elliptical corner crack in a finite plate in tension [Newman79]
eft-newman84         a/d, a/c, c/b, phi             elliptical crack in a finite plate in tension [Newman79]
sset-tada73          a/d                            single sided edge crack in a plate in tension [Tada73]
dset-tada73          a/d                            double sided edge crack in a plate in tension [Tada73]
compact-tada73       a/d, depth, width              compact specimen in tension (scale is in load units not stress units)  [Tada73]
ct-fedderson66       a/d                            centre cracked plate in tension [Fedderson66]
ct-koiter65          a/d                            centre cracked plate in tension [Koiter65]
qcct-mcdonald07      a/d                            vertically constrained coupon with corner crack in tension [McDonald07]
ssht-bowie56         a/r                            single sided through crack in a circular hole in tension [Bowie56]
dsht-bowie56         a/r                            double sided crack through in a circular hole in tension [Bowie56]
dccht-newman81       a/d, a/c, c/b, a/r, phi        double sided corner crack in a hole in tension [Newman81]
serbb-shin04         a/d, a/c                       semi-elliptical surface crack in a round bar in bending [shin04]
serbb-murakami87     a/d, a/c                       semi-elliptical surface crack in a round bar in bending [Murakami87]
serbt-murakami87     a/d, a/c                       semi-elliptical surface crack in a round bar in tension [Murakami87]
serbb-murakami86     a/d, a/c                       semi-elliptical surface crack in a round bar in bending [Murakami86]
serbt-murakami86     a/d, a/c                       semi-elliptical surface crack in a round bar in tension [Murakami86]
esb-murakami87       a/d                            edge crack in a strip in bending [Murakami87]
est-murakami87       a/d                            edge crack in a strip in tension [Murakami87]
file:FILE            a/d, a/c                       read FILE for beta values 

# Crack growth Models

1. rainflow  Crack growth is calculated from rainflow cycles i.e., the
          stress intensity range comes from the range of the rainflow
          cycles. Note this has a slight re-ordering effect that may upset
          the order of any image plots created.

2. tension   Crack growth calculated from tension part of cycle i.e. from a valley to the next peak.


# da/dN data

The da/dN model consists of EQUATION:material where the
equation variable specifies the name of the da/dN equation and is
one of  [nasgro, paris, forman, walker, burchill, hartman, white, file]
The material varable specifies the name of the parameters to use for
that equation. If the values are given in --parameters they will
be used instead of the standard library values.

Name                                Ref.                 Coefficients
paris:default                       [none]               1.00000e-10 3.00000e0 
walker:default                      [none]               1.00000e-10 5.00000e-1 3.00000e0 
forman:default                      [none]               1.00000e-10 3.00000e0 6.00000e1 
nasgro:default                      [nasgro:aa7050t7451-LT, NASGR04.0] 3.00000e-1 2.00000e0 3.51600e1 8.00000e-1 2.20000e0 1.00000e-1 1.00000e0 1.00000e0 6.35000e-10 2.50000e0 3.81000e-5 
nasgro:aa7050t7451-LT               [Forman05]           3.00000e-1 2.00000e0 3.51600e1 8.00000e-1 2.20000e0 1.00000e-1 1.00000e0 1.00000e0 6.35000e-10 2.50000e0 3.81000e-5 
burchill:default                    [none]               1.00000e-10 3.00000e0 1.00000e-10 3.00000e0 
kujawski:default                    [none]               1.00000e-10 3.00000e0 5.00000e-1 
hartman:default                     [none]               1.00000e-10 1.00000e0 3.00000e1 3.00000e0 
paris:newman-aa7050t7451            [none]               1.59300e-11 3.66800e0 
forman:aa2024t3-sheet               [Schwarmann86]       7.13000e-9 2.70000e0 7.13000e1 
forman:aa2024t351-plate             [Schwarmann86]       5.00000e-9 2.88000e0 6.32000e1 
forman:aa2024t4-sheet               [Schwarmann86]       8.57000e-9 2.60000e0 5.81000e1 
forman:aa2024t6-sheet               [Schwarmann86]       2.00000e-8 2.62000e0 6.98000e1 
forman:aa2024t8-sheet               [Schwarmann86]       1.33000e-8 2.65000e0 6.53000e1 
forman:aa2024t851-plate             [Schwarmann86]       7.72000e-9 2.78000e0 6.14000e1 
forman:aa2219t851-plate             [Schwarmann86]       4.84000e-8 2.16000e0 5.75000e1 
forman:aa2618t6-sheet               [Schwarmann86]       8.56000e-9 2.58000e0 4.59000e1 
forman:aa6061t6-sheet               [Schwarmann86]       2.27000e-7 6.01000e1 1.66000e0 
forman:aa6061t651-plate             [Schwarmann86]       9.60000e-8 1.84000e0 4.12000e1 
forman:aa7010t73651-plate           [Schwarmann86]       2.06000e-8 2.46000e0 4.60000e1 
forman:aa7050t7352-forging          [Schwarmann86]       2.75000e-9 3.29000e0 6.40000e1 
forman:aa7050t73651-plate           [Schwarmann86]       4.11000e-9 2.98000e0 5.50000e1 
forman:aa7075t6-sheet               [Schwarmann86]       1.37000e-8 3.02000e0 6.39000e1 
forman:aa7075t7351                  [Schwarmann86]       6.27000e-9 2.78000e0 5.58000e1 
forman:aa7175t3652-forging          [Schwarmann86]       2.61000e-9 2.91000e0 3.80000e1 
forman:aa7178t651-plate             [Schwarmann86]       3.74000e-8 2.06000e0 3.07000e1 
forman:aa7475t7351-plate            [Schwarmann86]       3.24000e-8 2.32000e0 7.82000e1 
forman:aa7475t76-sheet              [Schwarmann86]       6.54000e-8 2.18000e0 7.99000e1 
forman:aa7475t7651-plate            [Schwarmann86]       9.30000e-9 2.73000e0 6.31000e1 
forman:a357t6-sandcasting           [Schwarmann86]       2.19000e-9 2.94000e0 4.15000e1 
forman:a357t6-investmentcasting     [Schwarmann86]       6.65000e-9 2.40000e0 3.82000e1 
hartman:jones13-aa7050t7451         [jones13]            7.00000e-10 1.00000e-1 4.70000e1 2.00000e0 
white:barter14-aa7050t7451          [white15]            2.54819e-1 1.10247e0 4.35832e0 2.30859e1 3.42017e-2 4.71784e-1 3.15400e1 
white:chan16-aa7050t7451            []                   2.91862e-1 1.26351e0 3.55283e0 2.22432e1 3.92409e-2 5.55131e-1 4.14592e1 
file:FILE                                                Read FILE of tabular dadn data.

# File formats

## Crack file

The crack growth file is in the following format:

<line> <block> <a>
...

or 

<block> <a>
...

or

<a>
...
                   
Blank lines in the file indicate non-contiguous measurements. If
<line> or <block> are missing the program will assume the readings are
one block apart with each block measured at line 0. Use the same
format for the entire file. Where <line> represents the corresponding
line no. (starting at 0) of the sequence file, and <block> is the
no. of the block at that crack depth. Strictly speaking, the actual
block numbers are not used by easigro with only the difference between
the block numbers in contiguous measurements used. Easigro only
matches the average crack growth rate using:

   rate = (growth between measurements) / (no. of blocks between measurements).

## Optimise file

Each line in the optimise file is a list of easigro command lines
(without the 'easigro' command) with each line containing the easigro
options that will best reproduce the crack growth calculation for the
associated crack growth curve that it is trying to match. Note: Only
the material model specified by the main command line that invokes the
optimisation will be used for all crack predictions, since those will
be the parameters that are optimised. Any other material
specifications will be ignored.

The format of the optimisation file is:

<easigro option> ... --crack <FILE1>
<easigro option> ... --crack <FILE2>
...

## Beta file

All lines beginning with a # are treated as a comment and ignored. The
format of the beta file is

\# Comment describing the contents of the file
a/d beta 
...

## Dadn file

All lines beginning with a # are treated as a comment and ignored. The
format of the file is:

\# Comment describing the contents of the file
r1 r2 ....
dadn1 deltaK1_r1 deltaK1_r2 ....
dadn2 deltaK2_r1 deltaK2_r2 ....
...


# References

[Newman79] J. C. Newman , Jr. and I. S. Raju
           Analyses of surface cracks in finite plates under tension or bending loads
           NASA Technical Paper 1578
           December 1979

[Newman81]  J. C. Newman Jr. and I. S. Raju
           Stress intensity factor equations for cracks
           in three-dimensional finite bodies, 
           NASA Technical Memorandum 83299, 1981 p 1--49.

[Newman81]  J. C. Newman Jr. and I. S. Raju
           Stress-intensity factor equations for
           cracks in three-dimensional finite bodies
           subjected to tension and bending loads, 
           NASA Technical Memorandum 85739, 1984.

[Anderson05] T. L. Anderson
           Fracture Mechanics - Fundamentals and Applications
           Taylor and Francis 3rd Edition 2005

[Tada73] H. Tada, P. C. Paris and G. R. Irwin
           The Stress Analysis of Cracks Handbook
           1973

[Murakami87] Y. Murakami
          Stress Intensity Factors Handbook. Vol 2
          Pergamon Press, Oxford, , 1987

[Murakami87a] Yukitaka Murakami and Hideto Tsuru
          Stress Intensity factor equations for a semi-elliptical surface crack in a shaft under bending
          1986

[Schwarmann86] L. Schwarmann
          Material Data of High-Strength Aluminium Alloys for Durability Evaluation of Structures
          Aluminium-Verlag 1986
          Note: The data from this report has been converted from mm/cycle to m/cyclic by factoring cf by 1e3.

[Fedderson66] 
         Taken from Damage Tolerant Design handbook from AFGROW documentation.

[Kujawski01] Daniel Kujawski, 
          A fatigue crack driving force parameter with load ratio effects
          International Journal of Fatigue, Vol 23, S239-S246, 2001

[Walker70] K. Walker
          The effect of stress ratio during crack propagation and fatigue for {2024-T3} and {7075-T6} aluminum
          Effects of Environment and Complex Load History for Fatigue Life, 
          American Society for Testing and Materials,Special Technical Publication 462, 1970

[Jones12] Jones, R., Molent, L. & Walker, K.  
          Fatigue crack growth in a diverse
          range of materials, International Journal of Fatigue Vol. 40,pages 43--50, 2012

[Hartman70] A. Hartman and J. Schijve
          The effects of environment and load frequency on the
          crack propagation law for macro fatigue crack growth in aluminum alloys,
          Engineering Fracture Mechanics, Vol. 1(4), 1970

[Shin04] C.S. Shin and C. Q. CAI
          Experimental and finite element analyses on stress intensity
          factors of an elliptical surface crack in a circular shaft under
          tension and bending,
          International Journal of Fracture 129: 239–264, 2004.

[Forman05] R. G. Forman, V. Shivakumar, J. W. Cardinal , L. C. Williams and P. C. McKeighan 
                        Fatigue Crack Growth Database for Damage Tolerance Analysis, 
                        DOT/FAA/AR-05/15, 2005.

