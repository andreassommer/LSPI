# LSPI

**L**ocal **S**hape-**P**reserving (C2-smooth) **I**nterpolation

(c) 2026, Andreas Sommer <code@andreas-sommer.eu>  
(c) 2024, 2025, 2026  Julian Niederer <julian.niederer@iwr.uni-heidelberg.de>

Code based on prototype implementation by Julian Niederer.

Literature:  
J. Niederer, A. Sommer, H. G. Bock, E. Kostina, C. Fleck, A. Flexeder, M. Bitzer, C. Bertsch:
"A local monotonicity-preserving C2-smooth interpolation method for multivariate gridded data with applications in the automotive industry" (2026)
in review

## Usage

### Compilation

Compile using make:
```console
make all
```

or manually
```console
g++ -std=c++11 -Wall -O2 -I. -o test_1D_local_interpolation test_1D_local_interpolation.cpp
g++ -std=c++11 -Wall -O2 -I. -o test_2D_local_interpolation test_2D_local_interpolation.cpp
```

Note that changing the optimization flag "-O2" significantly influences runtime.

### Running

There are two exemplary test programs, for which help can be displayed via `--help` flag.

```console
./test_1D_local_interpolation --help
Usage:  test_1D_local_interpolation.exe n nrep output randseed
          n: number of evaluation points between grid points  [default if not specified: 10]
       nrep: number of repetitions for time measurement       [default if not specified: 1000000]
     output: 0: display  runtime  -- 1: display also values   [default if not specified: 0]
   randseed: random seed for shuffling evaluation points      [default value 0 disables randomization]
```

```console
./test_2D_local_interpolation --help
Usage:  C:\msys64\home\asommer\source\LSPI\test_2D_local_interpolation.exe silent n nrep
      n: number of evaluation points between grid points [default if not specified: 10]
   nrep: number of repetitions                           [default if not specified: 1]
 silent: 0: display only runtime --- 1: display interpolated values (Matlab format);
```

For further questions, contact the authors.