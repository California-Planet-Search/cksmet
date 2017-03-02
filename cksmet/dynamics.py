"""
Module for computing properties relating to the dynamics of different planets.
"""

def add_delta(pairs):
    """
    Compute the seperation between adjacent planet pairs interms of
    mutual hill-radii.

    """
    
    pairs['inner_mass'] = radius_to_mass(pairs.inner_iso_prad)
    pairs['outer_mass'] = radius_to_mass(pairs.outer_iso_prad)
    _delta = delta(
        pairs.inner_mass,
        pairs.outer_mass, 
        pairs.mass, 
        pairs.inner_smax, 
        pairs.outer_smax
    )
    pairs['delta'] = _delta
    return pairs

def hill_radius(mass1, mass2, stellar_mass, a1, a2):
    """
    mass1 (float): mass of inner planet (Earth Masses)
    mass2 (float): mass of outer planet (Earth Masses)
    stellar_mass (float): stellar mass (Solar Masses)
    a1 (float): semi-major axis (AU)


    """
    mass1 = np.array(mass1) * u.M_earth
    mass2 = np.array(mass2) * u.M_earth

    mass1 = mass1.to(u.M_sun).value
    mass2 = mass2.to(u.M_sun).value

    _hill_radius = ( 
        ((mass1 + mass2)/(3*stellar_mass))**(1.0/3.0) * 
        0.5*(a1 + a2)
    )
    return _hill_radius

def delta(mass1, mass2, stellar_mass, a1, a2):
    _hill_radius = hill_radius(mass1, mass2, stellar_mass, a1, a2)
    _delta = (a2 - a1)/_hill_radius
    return _delta

def radius_to_mass(radius):
    """
    Implement Weiss-Marcy Mass radius relationship
    """

    flux = 100 * FLUX_EARTH
    
    if isinstance(radius,Iterable):
        mass = map(radius_to_mass,radius)
        mass = np.array(mass)
        return mass
    
    if radius < 1.5: 
        rho = 2.43 + 3.39 * radius # g/cc WM14
        mass = (rho / 5.51) * radius**3
    elif 1.5 <= radius and radius <= 4.0: 
        mass = 2.69 * radius**0.93
    elif 4.0 < radius and radius < 9.0: 
        mass = 0.298 * radius**1.89 * flux**0.057
    elif 9.0 <= radius:
        mass = np.nan
    else:
        mass = np.nan
    return mass




