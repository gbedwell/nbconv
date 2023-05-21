# nbconv 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.

# nbconv 1.0.0

* nbconv is on CRAN!

# nbconv 1.0.1

* Fixed bug in nb_sum_saddlepoint() that could lead to NaN errors when solving the saddlepoint equation with uniroot.
* Added parallelization in rnbconv() to speed up large-scale sampling.
* Updated notation in nb_sum_saddlepoint() for clarity.