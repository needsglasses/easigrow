# Crack growth calculation

Here is an example of a simple crack growth calculation for a
compact tension coupon. The compact tension geometry factor uses the
applied load on the coupon rather than the far-field stress, so for
this example the scaling factor is 5000 N and because consistent units
are MN this becomes 5000E-6 MN. The forward distance of 50 mm
represents the distance from the centre-line of the pin attachment
holes to the free edge but the sideways distance is actually the full
width of the coupon of 12.5 mm. The initial starting distance of 15 mm
is the length of the crack from the centre-line of the pin
attachments. All distance units are given in m.

```
./easigrow -s 5000e-6 -b compact-tada73 -o block,a,k -a 0.015 -e 0.05  
          -q sequences/seq2.txt --forward=0.05 --sideways=0.0125 -n 10 

#  easigrow: version 0.0.0
#  
#  Options: 
#  a: [0.015]
#  a_limit: [0.05]
#  block_limit: 1000.0
#  params: []
#  output_vars: ["block", "a", "k"]
#  output_every: 10
#  output_lines: [1]
#  scale: 0.005
#  beta: "compact-tada73"
  component: Component { forward: 0.05, sideways: 0.0125, radius: inf, material: Properties { yield_stress: 450.0, k1c: 33.0, youngs_modulus: 71000.0 } }
#  seq_infile: "../check/sequences/seq1.txt"
#  seq_mods: SequenceModifiers { cap_max: None, cap_min: None, remove_bigger: None, remove_smaller: None, cycles: false, reorder: false, turning_points: false, outfile: None }
#  cycle_method: Rainflow
#  cycle_infile: ""
#  cycle_mods: CycleModifiers { cap_max: None, cap_min: None, remove_smaller: None, remove_bigger: None, remove_region: None, outfile: None }
#  dadn: "white:barter14-aa7050t7451"
#  No parameters given, obtaining from material library for white:barter14-aa7050t7451
#  da/dn equation: White { a: 0.25481858, b: 1.10247048, c: 4.35831677, d: 23.08586582, e: 0.03420171, f: 0.4717843, kic: 31.54, cite: "unknown", units: "m" }
#  da/dN (m) = exp[(2.5481858e-1 * ΔKeff^3 - 1.10247048e0 * ΔKeff^2 + 4.35831677e0 * ΔKeff - 2.308586582e1) 
#              + (dkic - ΔK)^-0.4717843] [White14:unknown]
#  where dkic = 3.154e1 (1 - R) and ΔKeff = ΔK / (1 - R)^0.03420171
       block            a            k 
      0.0000  1.500000e-2       0.0000 
     10.0000  1.543065e-2      10.2813 
     20.0000  1.589049e-2      10.5276 
     30.0000  1.638381e-2      10.7977 
     40.0000  1.691594e-2      11.0965 
     50.0000  1.749375e-2      11.4305 
     60.0000  1.812628e-2      11.8087 
     70.0000  1.882580e-2      12.2439 
     80.0000  1.960970e-2      12.7555 
     90.0000  2.050391e-2      13.3739 
    100.0000  2.155006e-2      14.1522 
    110.0000  2.282240e-2      15.1927 
    120.0000  2.447781e-2      16.7376 
    130.0000  2.697975e-2      19.6222 
    138.4711  1.744446e-1      60.7851 
#  Failure Event: a[0.17444463103533853] >= a_limit[0.05]
#  Failure Event: a[0.17444463103533853] > depth[0.05]
#  Failure Event: k[60.78514498903261] > k1c[33]
```
