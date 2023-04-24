# NB parameters

var( ' m r t' )
p = r / ( m + r )
q = m / ( m + r )
term1 = p * ( m + r )
x = 1 - ( ( m * exp(t) ) / ( m + r ) )
z = x * ( m + r )
term2 = z.full_simplify()

print("CGF term 1")
print(term1)
print("")

print("CGF term 1")
print(term2)
print("")

K = r * ( log(term1) - log(term2) )

assume( term2 > 0 )

dK = factor( derivative(K, t, 1).full_simplify() )
ddK = factor( derivative(K, t, 2).full_simplify() )

print("CGF")
print(K)
print("")

print("CGF first derivative")
print(dK)
print("")

print("CGF second derivative")
print(ddK)
print("")

reset()

var( ' p r t ' )

M(t) = ( ( p ) / ( 1 - ( 1 - p ) * exp(t) ) )^r

dM(t) = factor( derivative( M(t), t, 1 ).full_simplify() )

ddM(t) = factor( derivative( M(t), t, 2 ).full_simplify() )

dddM(t) = factor( derivative( M(t), t, 3 ).full_simplify() )

ddddM(t) = factor( derivative( M(t), t, 4 ).full_simplify() )

M1 = dM(t).subs(t == 0)

M2 = ddM(t).subs(t == 0)

M3 = dddM(t).subs(t == 0)

M4 = ddddM(t).subs(t == 0)

print("M1")
print(M1)
print("")

print("M2")
print(M2)
print("")

print("M3")
print(M3)
print("")

print("M4")
print(M4)
print("")



# cumulants of MGF

var( 'phi q p t mu' )

## MGF in terms of p
print( "CGF in terms of p and phi" )
print( "" )

M = -phi * log( 1 - ( ( 1 - p ) / p ) * ( exp(t) - 1 ) )
print( "CGF:" )
print( M )
print( "" )

print( "First derivative of CGF:" )
dM = factor( derivative( M, t, 1 ).full_simplify() )
print( dM )
print( "" )

print( "First cumulant:" )
k1 = dM.subs(t == 0)
print( k1 )
print( "" )

print( "Second derivative:" )
ddM = factor( derivative( M, t, 2 ).full_simplify() )
k2 = ddM.subs(t == 0)
dddM = factor( derivative( M, t, 3 ).full_simplify() )
k3 = dddM.subs(t == 0)
ddddM = factor( derivative( M, t, 4 ).full_simplify() )
k4 = ddddM.subs(t == 0)

## MGF in terms of mu
M = -phi * log( 1 - ( ( 1 - ( phi / ( phi + mu ) ) ) / ( phi / ( phi + mu ) ) ) * ( exp(t) - 1 ) )
dM = factor( derivative( M, t, 1 ).full_simplify() )
k1 = dM.subs(t == 0)
ddM = factor( derivative( M, t, 2 ).full_simplify() )
k2 = ddM.subs(t == 0)
dddM = factor( derivative( M, t, 3 ).full_simplify() )
k3 = dddM.subs(t == 0)
ddddM = factor( derivative( M, t, 4 ).full_simplify() )
k4 = ddddM.subs(t == 0)

k1
k2
k3
k4