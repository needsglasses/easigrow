The 'data/sequences' directory contains simple loading sequences which
can be used for testing easigrow. They have been previously used for
exploring fatigue crack growth effects which are described in the two
papers below.

In summary, the paper on rainflow counting showed that there was
essentially no memory effect. So while reordering a sequence using
rainflow counting works for crack growth, any process that reproduces
the same crack extension will result in the same growth, based on
crack extension occuring on the incresing loading part of the cycle
e.g. if we have 10 repeated cycles going from 0.0 to 1.0. We will get
the same growth if we have 10 cycles of 0--0.5 followed by 10 cycles
of 0.5--1.0 . It seems that the stress ratio effect is more signficant
for short periods of growth than is indicated from constant amplitude
testing at higher stress ratios. 

Secondly, the measurements on crack closure showed that the variations
in crack grwoth rate at different stress ratios did not correlate with
measurements of crack closure. Crack closure was very slow to build
up, occuring over a distance of millimeters, whereas the change in
growth rate due to different mean stress occured within a few cycles.

The papers are:

'''
@Article{white14,
  author = 	 {Paul White and David S. Mongru},
  title = 	 {Fractographic study on the use of rainflow counting for small and long
              cracks in {AA7050}},
  journal =  {Advanced Materials Research},
  year = 	 2014,
  volume = 	 {891--892},
  pages = 	 {687--692}
}

@Article{white18,
  author = 	 {P. D. White and S. A. Barter and N. Medhekar},
  title = 	 {Comparison of fractographic and strain gauge  measurement of
              load ratio effects under variable amplitude loading},
  journal =  {International Journal of Fatigue},
  volume =   112,
  pages =    {240--252},
  month =    {July},
  year = 	 2018
}
'''
