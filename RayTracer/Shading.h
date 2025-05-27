#pragma once

/*
surface shading:
at each surface point we just have 1 sample from each light source

right now just a basic lambert surface: the brdf is (r/pi)

so in the integral: (r/pi)*lightcolor*cos(theta) 

R can just be total refletance percentage, cos(theta) is angle between normal and light, 4pi comes from surface area of halfsphere



ray pathing:
every time it hits an object it will create a reflection raw (perfect reflection), refract (snells law) and their percentages will be 
based on fresnell giving percentage of reflection/refraction. 
Then there will be 2-3 bounces


*/