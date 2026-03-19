### Radiometry

Radiometry provides a set of ideas and mathematical tools to describe light propagation and reflection. It forms the basis of the derivation of the rendering algorithms that will be used throughout the raytracer. 

> Interestingly enough, radiometry was not originally derived from first principles using the physics of light but was built on an abstraction of light based on particles flowing through space. As such, effects like polarization of light do not naturally fit into this framework, although connections have since been made between radiometry and Maxwell’s equations, giving radiometry a solid basis in physics

*Radiative transfer* is the transfer of radiant energy (duh). It is based on radiometric principles and operates at the _geometric optics_ level, where macroscopic properties of light suffice to describe how light interacts with objects much larger than the light’s wavelength. We assu

We assume the following:
- _Geometric optics model_ - Going past geometric optics includes sampling objects the same size of the wavelength, or even finer like quantum mechanics to describe the interactions with atoms. Instead, we operate on the macroscopic properties of light suffice to describe how light interacts with objects much larger than the light’s wavelength.
- _Linearity_ - The combined effect of two inputs to an optical system is always equal to the sum of the effects of each of the inputs individually.
- _Energy conservation_ - When light scatters from a surface or from participating media, the scattering events can never produce more energy than they started with.
- _No polarization_ - Electromagnetic radiation including visible light is _polarized_, but curiously, this additional polarization state of light is essentially imperceptible to humans without additional aids like specialized cameras or polarizing sunglasses. Therefore, the only relevant property of light is its distribution by wavelength.
- _No fluorescence or phosphorescence_ - The behavior of light at one wavelength is completely independent of light’s behavior at other wavelengths or times. 
- _Steady state_ - Light in the environment is assumed to have reached equilibrium, so its radiance distribution is not changing over time.

### Basic Quantities
There are four radiometric quantities that are central to rendering: [[Energy]], [[Flux]], [[Irradiance and Radiant Exitance| irradiance / radiant exitance]], [[Intensity]], and most importantly, [[Radiance]]. They can each be derived from energy by successively taking limits over time, area, and directions. All of these radiometric quantities are in general wavelength dependent.

### Luminance and Photometry
_Photometry_ is the study of visible electromagnetic radiation in terms of its perception by the human visual system.

All the radiometric measurements like flux, radiance, and so forth have corresponding photometric measurements. Each spectral radiometric quantity can be converted to its corresponding photometric quantity by integrating against the spectral response curve $V(\lambda)$, which describes the relative sensitivity of the human eye to various wavelengths. 

_Luminance_ measures how bright a spectral power distribution appears to a human observer. 

We will denote luminance by $Y$; it is related to spectral radiance by 
$$Y=\int_\lambda L(\lambda)V(\lambda)d\lambda$$

Luminance and the spectral response curve $V(\lambda)$ are closely related to the [[XYZ]] representation of color. The units of luminance are candelas per meter squared.