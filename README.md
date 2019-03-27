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
- Generates a psuedo image of the crack growth patterns

## Getting started

Easigrow provides both a command line tool for running simulations;
and a Rust library for writing your own algorithms using the
pre-existing building blocks.

Here is a simple example of growing a crack from the command line:

```bash
./target/release/easigrow -q sequences/seq1.txt -s 300 -n 100

#  easigrow: version 0.0.0
#  
#  Options: 
#  a: [0.00001, 0.00001]
...
       block            a 
      0.0000  1.000000e-5 
    100.0000  1.514032e-5 
    200.0000  4.002745e-5 
    300.0000  6.635652e-4 
    301.5985  7.468295e-4 
#  Failure Event: a[0.0007468294682320596, 0.001000259809490663] >= a_limit[0.001, 0.001]
```

### Documentation

There are two sources of documentation:

- Run `easigrow --help`
  if you want to use Easigrow as a command line tool
  without writing code.
- Run `Cargo doc` to generate documentation for all the library public functions and use Easigrow as a library inside your own code.

### Installation as a command line tool

You will need a stable Rust compiler, [grab one][Rust] if you do not have one
yet. Then, you can download the code, build it and install it by running:

```bash
cargo install --git https://github.com/needsglasses/easigrow
```

This will produce the a `easigrow` binary in `~/.cargo/bin`.

### Usage as a library

You can add Easigrow as a dependency in your project's `Cargo.toml`:

```toml
[dependencies]
easigrow = {git = "https://github.com/needsglasses/easigrow"}
```

## Contributing

If you want to contribute to Easigrow, there are several ways to go:
improving the documentation; testing the code on your systems to find
bugs; adding new algorithms or data; providing feature requests.

See the [AUTHORS](AUTHORS) file for a list of contributors to the code.

## License

This software is licensed under the MIT license, see the
[LICENSE](LICENSE) file for legal text.

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in the work by you, shall be licensed under the same MIT license,
without any additional terms or conditions.

[Rust]: https://www.rust-lang.org/downloads.html
[issues]: https://github.com/needsglasses/easigrow/issues/new
[contributing]: Contributing.md
[user_manual]: http://needsglasses.github.io/easigrow/latest/book/
[devdoc]: http://needsglasses.github.io/easigrow/latest/easigrow/
