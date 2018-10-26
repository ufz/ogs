
import OpenGeoSys
from math import pi, sin

a = 2.0*pi
b = 2.0*pi

def solution(x, y):
    return sin(a*x-pi/2.0) * sin(b*y-pi/2.0)

# - laplace(solution) = source term
def laplace_solution(x, y):
    return a*a * sin(a*x-pi/2.0) * sin(b*y-pi/2.0) + b*b * sin(a*x-pi/2.0) * sin(b*y-pi/2.0)

# source term for the benchmark
class SinXSinYSourceTerm(OpenGeoSys.SourceTerm):
    def getFlux(self, t, coords, primary_vars):
        x, y, z = coords
        value = laplace_solution(x,y)
        Jac = [ 0.0 ]
        return (value, Jac)

# instantiate source term object referenced in OpenGeoSys' prj file
sinx_siny_source_term = SinXSinYSourceTerm()
