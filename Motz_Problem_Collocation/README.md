# Motz_Problem_Collocation

Octave code for obtaining coefficient of the Motz problem using collocation method.

My next goal is to implement the Python code using either numpy.trapz or scipy.integrate and one of its variants (quad, fixed_quad (guass quadrature), romb (romberg), simps (composite simpson's rule)). 

Another thing would be to use tensorflow to the motz problem for a fast publication? possibly? The dataset would be (x,y,z) triples. I spent some time learning dense networks implementation & although I don't have a GPU, I managed to run it on my CPU. So it would be a good thing to try this at least. 

Another possibility is to use exact integration simplifications & using integration by parts to simplify the process though this may be very difficult. 

I don't even know if this problem is worth my time. This problem is relatively dead & I don't want to spend time on enriched finite difference methods. 

But I will at least try to do this in python. I will compare integration in scipy with manual integration using numpy+numba & see which one is faster. 

Another possibility is to use crude mesh to get a first estimate of the coefficients, then use the coefficients in the formula itself & see if this scaling helps things? So for example I can start with 5 coefficients, calculate their values (aj). In the next step, I will first predict a(j+1) by extrapolation. Then use aj.bj as coefficients, wherein bj is the unknown correction (in the order of 1) and aj is the either known value from previous step or the extrapolated one. 
