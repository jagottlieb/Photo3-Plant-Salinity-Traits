"""Maximum stem water storage per unit leaf area (VWT) for the Pecan lab trees."""

import math

STEM_HEIGHT = 1.5           # stem height (m); average lab tree, treated as a cylinder (no taper)
STEM_DIAMETER = 0.01        # stem diameter (m); average lab tree
WOOD_DENSITY = 0.60         # pecan wood density (g/cm3); green-volume specific gravity, USFS (2009)
CELL_WALL_DENSITY = 1.54    # density of wood cell wall material (g/cm3); Siau (1984), Christoffersen et al. (2016)
LEAF_AREA = 0.36            # total leaf area (m2); average lab tree (3600 cm2)


def vwt(height: float, diameter: float, wood_density: float, cell_wall_density: float, leaf_area: float) -> float:
    """Saturated stem water volume per unit leaf area (m3 water / m2 leaf).

    Assumes the whole stem cross-section is sapwood and its pore space (1 - WD / cell wall density)
    is water-filled at saturation.
    """
    stem_volume = math.pi * (diameter / 2) ** 2 * height
    water_fraction = 1 - wood_density / cell_wall_density
    return stem_volume * water_fraction / leaf_area


print(f"VWT = {vwt(STEM_HEIGHT, STEM_DIAMETER, WOOD_DENSITY, CELL_WALL_DENSITY, LEAF_AREA):.6f} m3/m2 leaf")
