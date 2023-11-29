# MARS calculator

# Python package

This Python package is a MARS calculator, based on Ferry et al.'s work [1].
It depends on:

* `islpy`,
* `numpy`

You need to have the `barvinok` software installed and `iscc` to be available
in your PATH. Get `barvinok` [here](https://barvinok.sf.net/).

The code is provided under the MIT license.

# Examples

The `examples` directory contains Python scripts computing the MARS for 
all the benchmarks found in the 
[MARS paper](https://impact-workshop.org/impact2023/papers/paper1.pdf). 

# IMPACT23 calculator

The calculator created for the IMPACT'23 submission and its definition of
MARS (not parametrized by the tile) is available in the `impact23` directory.
MARS computed using the `mars` package are parametrized by the tile they belong
to.

# References

\[1\] Corentin Ferry, Steven Derrien, and Sanjay Rajopadhye. « Maximal Atomic 
irRedundant Sets: a Usage-based Dataflow Partitioning Algorithm ». In: 13th 
International Workshop on Polyhedral Compilation Techniques (IMPACT’23). 2023. 
[Paper](https://impact-workshop.org/impact2023/papers/paper1.pdf) 
[Slides](https://impact-workshop.org/impact2023/slides/slides1.pdf)
