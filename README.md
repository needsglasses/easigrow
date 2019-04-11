# Easigrow
Create models of fatigue crack growth

Easigrow is a fatigue crack growth program that specialises in
calculating model coefficients that best fit crack growth data.

Easigrow is based on the Australian Defence Science and Technology report:
A guide to the program easigro for generating optimised fatigue crack
growth models, Paul White, DST-Group-TR-3566, Feb. 2019.

Note the addition of 'w' here for the code on Github.
The version number has also been reset to 0.0.0. 

## Features

- Calculates the crack growth size
- Choice of beta models, material properties and crack growth equations
- Optimises model coefficients to match data
- Generates a pseudo image of the crack growth patterns

## Getting started

Easigrow provides both a command line tool for running simulations
and a Rust library for writing your own algorithms using the
pre-existing building blocks.

Here is a simple example of growing a crack from the command line:

```
easigrow -q data/rainflow-seq2.txt -s 300 -n 100

#  easigrow: version 0.0.0
#  
#  Options: 
#  a: [0.00001, 0.00001]
...
       block            a 
      0.0000  1.000000e-5 
    100.0000  2.189452e-5 
    200.0000  5.234739e-4 
    204.0975  7.473615e-4 
#  Failure Event: a[0.0007473615299278992, 0.0010002610857908832] >= a_limit[0.001, 0.001]
```

### Documentation

In addition to the above mentioned report, there are the following
sources of documentation:

- Run `easigrow --help`
  if you want to use Easigrow as a command line tool
  without writing code.
- Examples using the command line in the `doc` directory
- Run `Cargo doc` to generate documentation for all the library public
  functions and use Easigrow as a library inside your own code.

### Installation as a command line tool

You will need a stable Rust compiler, [grab one][Rust] if you do not have one
yet. Then, you can download the code, build it and install it by running:

```bash
cargo install --git https://github.com/needsglasses/easigrow
```

This will produce the `easigrow` binary in `~/.cargo/bin`.

### Usage as a library

You can add Easigrow as a dependency in your project's `Cargo.toml`:

```toml
[dependencies]
easigrow = {git = "https://github.com/needsglasses/easigrow"}
```

## Contributing

Contributions to Easigrow are welcome. There are several ways to help:
improving the documentation; testing the code on your systems to find
bugs; adding new algorithms or data; providing feature requests.

See the [AUTHORS](AUTHORS) file for a list of contributors to the code.

## License

This software is licensed under the MIT license, see the
[LICENSE-MIT](LICENSE-MIT) file for legal text.

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, shall be licensed under the same MIT license,
without any additional terms or conditions.

[Rust]: https://www.rust-lang.org/downloads.html
[issues]: https://github.com/needsglasses/easigrow/issues/new
[user_manual]: http://needsglasses.github.io/easigrow/latest/book/
[devdoc]: http://needsglasses.github.io/easigrow/latest/easigrow/

## Other Programs

Fatigue crack growth programs

* [Afgrow] USAF crack growth software for Windows
* [Nasgro] NASA crack growth software with large material library
* [Fastran] Crack growth using closure in Fortran
* [Kujawski] Interactive web based software for crack initiation and growth

and some others on Github that may need further work:

* [pdprop] Simulation of metal fatigue crack propagation with accounting for material memory effects in Fortran
* [FCGR] Fatigue Crack Growth in Reliability using R
* [metal-fatigue] Start of a fatigue library in Python with some rainflow counting routines

[Kujawski]: https://scholarworks.wmich.edu/cgi/viewcontent.cgi?article=1057&context=fracaa
[efatigue]: https://www.efatigue.com/
[Fastran]: https://www.researchgate.net/publication/24314111_FASTRAN-2_A_fatigue_crack_growth_structural_analysis_program
[Nasgro]: https://www.swri.org/consortia/nasgro
[Afgrow]: https://www.afgrow.net/
[metal-fatigue]: https://github.com/alexm2303/metal-fatigue
[pdprop]: https://github.com/pdprop/pdprop
[FCGR]: https://github.com/cran/FCGR
