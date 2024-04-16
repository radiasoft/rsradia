# rshyst
The `rshyst` package provides tools for modeling hysteresis dynamics in magnetic materials.

## Installation
The `rshyst` package can be installed by running the command:

```
python3 -m pip install ./rshyst
```

from one directory above the downloaded repository and the *pyproject.toml* file.

## Documentation

The `rshyst` documentation, still under development, can be found at:

* https://github.com/radiasoft/rsradia/wiki/rshyst

## References

The models and methods used in the `rshyst` package are based on 

### Jiles-Atherton Model

The Jiles-Atherton model implemented in the package is based on the equations first presented in:

* D. Jiles, D. Atherton, "Theory of ferromagnetic hysteresis," *Journal of applied physics* **55** (1984)
  * [DOI: 10.1063/1.333582](https://doi.org/10.1063/1.333582)

with the extensions for anisotropic dynamics given later in:

* A. Ramesh, D. Jiles, J. Roderick, "A model of anisotropic anhysteretic magnetization," *IEEE Transactions on Magnetics* **32** (1996)
  * [DOI: 10.1109/20.539344](https://doi.org/10.1109/20.539344)

A number of the numerical implementations used in the package (particularly the numerical integration schemes) were chosen based on results and recommendations presented in:

* R. Szewczyk, "Computational problems connected with Jiles-Atherton model of magnetic hysteresis," *Recent Advances in Automation, Robotics and Measuring Techniques*, Springer International Publishing (2014)
  * [DOI: 10.1007/978-3-319-05353-0_27](https://doi.org/10.1007/978-3-319-05353-0_27)

### Preisach Model

The Preisach model implemented in this package is based on the original proposed by Preisach in:

* F. Preisach, "Über die magnetische Nachwirkung," *Zeitschrift für physik* **94** (1935)
  * [DOI: 10.1007/BF01349418](https://doi.org/10.1007/BF01349418)

## License

Copyright 2024 RadiaSoft LLC

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
