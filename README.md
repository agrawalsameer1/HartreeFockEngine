# HartreeFockEngine
A C++ implementation of the Hartree Fock algorithm for simulating quantum n-body systems, using the STO-6G basis set.

Currently, algorithm can accurately calculate ground-state energy and atomic radius for Hydrogen, Helium, and Beryllium.

**Calculated ground state energy values vs ground truth:**
|    | **Calculated Value** | **Ground Truth Value** |
| ------------- | ------------- | ------------- |
| Hydrogen | -0.499 Hartrees | -0.5 Hartrees |
| Helium | −2.860 Hartrees | -2.903 Hartrees |
| Beryllium | -14.566 Hartrees | -14.698 Hartrees |


**Calculated atomic radius values vs ground truth:**
|    | **Calculated Value** | **Ground Truth Value** |
| ------------- | ------------- | ------------- |
| Hydrogen | 52.8 picometers | 53 picometers |
| Helium | 30.5 picometers | 31 picometers |
| Beryllium | 128.9 picometers | 112 picometers |

Algorithm can also accurately calculate ground-state energy values for various ions:

|    | **Calculated Value** | **Ground Truth Value** |
| ------------- | ------------- | ------------- |
| H<sup>1-</sup> | -0.470 Hartrees | -0.519 Hartrees |
| He<sup>1+</sup> | −1.998 Hartrees | -1.998 Hartrees |
| Li<sup>2+</sup> | −4.480 Hartrees | -4.455 Hartrees |
| Li<sup>1+</sup> | −7.223 Hartrees | -7.268 Hartrees |
| Li<sup>1-</sup> | -7.381 Hartrees | -7.456 Hartrees |
| Be<sup>3+</sup> | -7.921 Hartrees | -7.997 Hartrees |
| Be<sup>2+</sup> | -13.479 Hartrees | -13.686 Hartrees |
