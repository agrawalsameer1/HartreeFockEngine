# HartreeFockEngine
A C++ implementation of the Hartree Fock algorithm for simulating quantum n-body systems, using the [STO-6G basis set]([url](https://www.basissetexchange.org/basis/sto-6g/format/nwchem/?version=1&elements=1,2,4)).

Currently, algorithm can accurately calculate ground-state energy and atomic radius for Hydrogen, Helium, and Beryllium.

**Calculated ground state energy values vs ground truth:**
|    | **Calculated Value** | **Ground Truth Value** | **Percent Error** |
| ------------- | ------------- | ------------- | ------------- |
| Hydrogen | -0.499 Hartrees | -0.5 Hartrees | 0.20% |
| Helium | −2.860 Hartrees | -2.903 Hartrees | 1.48% |
| Beryllium | -14.566 Hartrees | -14.698 Hartrees | 0.90% |


**Calculated atomic radius values vs ground truth:**
|    | **Calculated Value** | **Ground Truth Value** | **Percent Error** |
| ------------- | ------------- | ------------- | ------------- |
| Hydrogen | 52.8 picometers | 53 picometers | 0.38% |
| Helium | 30.5 picometers | 31 picometers | 1.61% |
| Beryllium | 128.9 picometers | 112 picometers | 15.1% |

Algorithm can also accurately calculate ground-state energy values for various ions:

|    | **Calculated Value** | **Ground Truth Value** | **Percent Error** |
| ------------- | ------------- | ------------- | ------------- |
| H<sup>1-</sup> | -0.470 Hartrees | -0.519 Hartrees | 9.44% |
| He<sup>1+</sup> | −1.998 Hartrees | -1.998 Hartrees | 0.0% |
| Li<sup>2+</sup> | −4.480 Hartrees | -4.455 Hartrees | 1.54% |
| Li<sup>1+</sup> | −7.223 Hartrees | -7.268 Hartrees | 0.62% |
| Li<sup>1-</sup> | -7.381 Hartrees | -7.456 Hartrees | 1.01% |
| Be<sup>3+</sup> | -7.921 Hartrees | -7.997 Hartrees | 0.95% |
| Be<sup>2+</sup> | -13.479 Hartrees | -13.686 Hartrees | 1.51% |
