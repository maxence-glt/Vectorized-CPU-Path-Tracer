

Light intensity is the power of light per unit area, representing the energy concentration of a light beam, often measured in 

- watts per square meter ($W/m^2)$
- lux ($lx$)
- lumens per square foot ($fc$)

It relates to brightness, generally decreasing with distance according to the inverse square law.

Consider now an infinitesimal light source emitting photons. If we center this light source within the unit sphere, we can compute the angular density of emitted power. _Intensity_, denoted by $I$, is this quantity. Over the entire sphere of directions, we have
$$I=\frac{Phi}{4\pi}$$
And for a differential cone of directions:
$$I=lim_{\Delta\omega\rightarrow0}\frac{\Delta\Phi}{\Delta\omega}=\frac{d\Phi}{d\omega}$$
To go back to power:
$$\Phi=\int_\Omega I(\omega)d\omega$$