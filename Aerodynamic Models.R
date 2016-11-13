#########################################################
### Functions for Estimating Avian Flight Performance ###
#########################################################

# These functions use morphometric data to estimate parameters related to flight performance such as the mechanical power required for flight based on aerodynamic models of flapping flight. The models are based on Pennycuick's "Modeling the flying bird" (2008, Accademic Press) with slight modifications explained in Claramunt & Wright (2016). Basic morphometric data required are mass, wingpann, and total wing area. Physical parameters such as gravity and air density have conventional default values but can be changed.

# References #

#Claramunt, S., & N. A. Wright. in press. Using museum specimens to study flight and dispersal. in M. S. Webster (ed.) Emerging Frontiers in Collections-based Ornithological Research: The Extended Specimen. Studies in Avian Biology.
# Pennycuick, C. J. 2008. Modeling the flying bird. Academic Press, Cambridge, MA.


### Arguments ###

# m		total body mass in kilograms
# B		wingspan in meters.
# Awing	total wing area in squared meters.
# V		velocity (m/s), true airspeed. If not provided, defaults to the minimum power velocity estimated from the power functions.
# Abody	frontal area of the body (m^2). If not provided, it is estimated from the body mass using the allometric equation 0.01*m^(2/3).
# Cbody	the body’s drag coefficient.
# Cpro	profile power constant
# rho	air density.
# g		gravitational acceleration.


### Minimum Power Velocity ###
# Is the velocity that minimizes the power required for flight.

V.mp <- function(m, B, Abody=0.01*m^(2/3), Cbody=0.1, rho=1.23, g=9.80665)
{(0.807*sqrt(m*g))/(sqrt(B*rho)*(Abody*Cbody)^0.25)}


### Induced Power ####
# Power required to support the weight of the bird

P.induced <- function(m, B, V, rho=1.23, g=9.80665) {(2*(m*g)^2)/(V*pi*rho*B^2)}

### Parasite Power ###
# Power required for counteracting the drag of the body.

P.parasite <- function(m, Abody=0.01*m^(2/3), Cbody=0.1, V, rho=1.23) {(Abody*Cbody*rho*V^3)/2}

### Profile Power ###
# The power required to overcome the drag of the wings that is not accounted for already as induced drag associated with lift.

P.profile <- function(m, Awing, B, Abody=0.01*m^(2/3), Cbody=0.1, Cpro=8.4, rho=1.23, g=9.80665) {
	
	Vmp <- V.mp(m=m, B=B, Cbody=Cbody, rho=rho, g=g)
	
	Pind <- P.induced(m=m, B=B, V=Vmp)
	Ppar <- P.parasite(m=m, Abody=Abody, V=Vmp, Cbody=Cbody, rho=rho)
	
	Ppro <- (Pind+Ppar)*Cpro*Awing/B^2
	
	return(Ppro)
}


### Total Mechanicas Power ###
# The sum of Induces, Parasite, and Profile powers.


P.mec <- function(m, Awing, B, Abody=0.01*m^(2/3), V=NULL, Cbody=0.1, Cpro=8.4, rho=1.23, g=9.80665) {
	
	if(is.null(V)) V <- V.mp(m=m, B=B, Cbody=Cbody, rho=rho, g=g)

	Pind <- P.induced(m=m, B=B, V=V)
	Ppar <- P.parasite(m=m, Abody=Abody, V=V, Cbody=Cbody, rho=rho)
	Ppro <- P.profile(m=m, Awing=Awing, B=B, Abody=Abody, Cbody=Cbody, Cpro=Cpro, rho=rho, g=g)
	return(Pind+Ppar+Ppro)
}

### Lift to Drag Ratio ###

Lift.Drag <- function(m, Awing, B, Abody=0.01*m^(2/3), V=NULL, Cbody=0.1, Cpro=8.4, rho=1.23, g=9.80665) {
	
	if(is.null(V)) V <- V.mp(m=m, B=B, Cbody=Cbody, rho=rho, g=g)

	Pind <- P.induced(m=m, B=B, V=V)
	Ppar <- P.parasite(m=m, Abody=Abody, V=V, Cbody=Cbody, rho=rho)
	Ppro <- P.profile(m=m, Awing=Awing, B=B, Abody=Abody, Cbody=Cbody, Cpro=Cpro, rho=rho, g=g)
		
	LDR <- m*g*V/(Pind+Ppar+Ppro)
	
	return(LDR)
}



### EXAMPLES ###

# Diffrenreces in flight efficiency among similar sized passerines

#Lift.Drag(m=0.017, B=0.324, Awing=0.014) # Barn swallow		

#Lift.Drag(m=0.0173, B=0.187, Awing=0.0087) # Rufous spinetail
	

# Mechanical power as a function of air speed for Anas penelope (Compare with figure 3.5 in Pennycuick, 2008).

#plot(8:25, P.mec(m=0.770, Awing=0.0829, B=0.822, V=c(8:25)), type="l", xlim=c(0,25), ylim=c(0,15), yaxs="i", xaxs="i", las=1, xlab="Airspeed m/s", ylab="Power W", lwd=2)

#lines(8:25, P.induced(m=0.770, B=0.822, V=8:25), col="blue3", lwd=2)
#lines(8:25, P.parasite(m=0.770, V=8:25), col="red3", lwd=2)
#lines(c(8,25), rep(P.profile(m=0.770, Awing=0.0829, B=0.822),2), col="green4", lwd=2)
#abline(v=V.mp(m= 0.770, B= 0.822), lty="dashed")

