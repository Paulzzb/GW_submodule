# How to use GW module

- 'git clone' the code from the repository to the same directory with the main program, or download the code from the repository and unzip it.
- run the demo_test to make sure the module is working properly.
- use the examples in the ./examples/ folder to start your own research.

Note: 
- The GW module is still under development, and some functions may not be available yet.
- The GW module ... 

# Structure of GW module
```
├── doc: documentations for GW module 
├── doc: documentations for GW module 
|   The following are classes in GW module.
|   ├── options: options for GW module, @GWOptions
|   The following are widely used functions.
|   ├── do_FFT & get_from_fftbox & put_into_fftbox:
|       a set of functions to do FFT.
|   ├── gvec_to_fft_index: N x N x N --> N^3, in kssolv manner.
├── doc: documentations for GW module 
├── src: source code for GW module, see ./src/README.md for detail
├── test: test functions for GW module.
├── example: simple examples for GW module
└── util: utility functions for GW module, mainly for format conversion and I/O conversion.
    see ./util/README.md for detail


```
hello