# Getting information about the loading sequence
Here is an example of the `--summary` option used to obtain
information about the sequence `seq2.txt` and the cycles it
contains. 

```
./easigro -q examples/seq2.txt --summary 

#  easiGro: version 0.0.0
#  

Sequence Summary
----------------

Source: examples/seq2.txt
Sequence mods: SequenceModifiers {
    cap_max: None,
    cap_min: None,
    remove_bigger: None,
    remove_smaller: None,
    cycles: false,
    reorder: false,
    turning_points: false,
    outfile: None
}
Length: 1340
Maximum: 1.0000e0 at line 3
Number of times maximum occurs: 160
Minimum: 0.0000e0 at line 0
Number of times minimum occurs: 160
Mean: 5.0000e-1
Number of points >= 0: 1340
Number of points < 0: 0
Number of non turning points: 0
Sequence: [0.0 0.9 0.1 1.0 0.0 ... 0.75 0.25 0.75 0.25 0.75]


Cycle Summary
-------------

Source: (Obtained from sequence using 'Rainflow' method)
Cycle mods: CycleModifiers {
    cap_max: None,
    cap_min: None,
    remove_smaller: None,
    remove_bigger: None,
    remove_region: None,
    outfile: None
}
Number of closed cycles: 667
Number of unclosed turning points: 6
Largest range: 1.0000e0, with valley of 0 at line 4, peak of 1 at line 3
Smallest range: 5.0000e-1 with valley of 0.25 at line 160, peak of 0.75 at line 161
Mean range: 6.9070e-1
Maximum R: 3.3333e-1
Minimum R: 0.0000e0
Mean R: 2.0025e-1
```
