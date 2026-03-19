Radiance is the radiant [[Flux|flux]] emitted, reflected, transmitted or received by a given surface, per unit solid angle per unit projected area. [[Irradiance and Radiant Exitance]] give us differential power per differential area at a point $p$ but not the directional distribution of power.

Radiance is defined by
$$\large L(p,\omega)=lim_{\Delta\omega\rightarrow0}\frac{\Delta E_\omega(p)}{\Delta\omega}=\frac{dE_\omega(p)}{d\omega}$$
where $E_\omega$ denotes irradiance at the surface that is perpendicular to the direction $\omega$. Radiance is not measured with respect to the irradiance incident at the surface $p$ lies on. In effect, this change of measurement area serves to eliminate the $cos\theta$ factor from [[Lambert's Cosine Law|Lambert’s law]] in the definition of radiance.

Radiance is the [[Flux|flux]] density per unit area, per unit solid angle. In terms of flux:
$$\large L=\frac{d^2\Phi}{d\omega A^\perp}$$
where $dA^\perp$ is the projected area of $dA$ on a hypothetical surface perpendicular to $\omega$. 
![[Pasted image 20260319122010.png]]

If radiance is given, then all the other values can be computed in terms of integrals of radiance over areas and directions. Another nice property of radiance is that it remains constant along rays through empty space. It is thus a natural quantity to compute with ray tracing.