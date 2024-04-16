"""constants.py

LICENSE STATEMENT

Copyright 2024 RadiaSoft LLC

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

MU0 = 1.25663706212e-6

UNITS = {"SI":1, "MM":1e-3, "TESLA": MU0, "GAUSS": MU0*1e4}

def convert(value, fromUnits, toUnits):
    """Converts the physical units of a magnetic quantity

    Args:
      * value: Quantitiy to be converted
      * fromUnits: Original units (str)
      * toUnits: Desired units (str)
    
    Note: Units must be one of 'SI', 'MM', 'TESLA', or 'GAUSS'
    """

    uFrom = UNITS[fromUnits.upper()]
    uTo  = UNITS[toUnits.upper()]

    return value * uFrom / uTo
