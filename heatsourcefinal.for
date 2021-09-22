	SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,JLTYP,TEMP,PRESS,SNAME)
	
	INCLUDE 'ABA_PARAM.INC' 
	
	DIMENSION FLUX(2), TIME(2), COORDS(3) 
	CHARACTER*80 SNAME 
	
	INTEGER kstep
	REAL t,x,y,z,pi
	REAL Q,a,b,bf,br,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,v,m,n,c
	REAL x7,y7,z7,x8,y8,z8,x9,y9,z9,x10,y10,z10,x11,y11,z11,x12,y12,z12,x13,y13,z13
	REAL x14,y14,z14,x15,y15,z15,x16,y16,z16

	   
		t = TIME(1)
		pi=3.414
		FLUX(1)=0
		FLUX(2)=0
!  speed of welding in y direction
		v=0.0067
!  coordinate
		x = coords(1)
		y = coords(2)
		z = coords(3)
!  welding arc 
		a = 0.005
!		b = 0.006
		c = 0.007
		bf = 0.005
		br = 0.020
!  for bead1 (transverse)		
		x1 = 0.495
		y1 = 0+v*t
		z1 = 0
		
		x2=x-x1
		y2=y-y1
		z2=z-z1
!  for bead2	(transverse)	
		x3 = 0.505
		y3 = 1.2-v*t
		z3 = 0
		
		x4 = x-x3
		y4 = y-y3
		z4 = z-z3
		
!  for bead3	(transverse)	
		x5 = 1.495
		y5 = 0+v*t
		z5 = 0
		
		x6 = x-x5
		y6 = y-y5
		z6 = z-z5		

!  for bead4	(transverse)	
		x7 = 1.505
		y7 = 1.2-v*t
		z7 = 0
		
		x8 = x-x7
		y8 = y-y7
		z8 = z-z7		
		
!  for bead5	(horizontal , heat sorce move in x-axis)	
		x9 = 2-v*t
		y9 = 0.295
		z9 = 0
		
		x10 = x-x9
		y10 = y-y9
		z10 = z-z9		
		
!  for bead6	horizontal , heat sorce move in x-axis)	
		x11 = 0+v*t
		y11 = 0.305
		z11 = 0
		
		x12 = x-x11
		y12 = y-y11
		z12 = z-z11		
		
!  for bead7	horizontal , heat sorce move in x-axis)	
		x13 = 2-v*t
		y13 = 0.895
		z13 = 0
		
		x14 = x-x13
		y14 = y-y13
		z14 = z-z13

!  for bead8	horizontal , heat sorce move in x-axis)	
		x15 = 0+v*t
		y15 = 0.905
		z15 = 0
		
		x16 = x-x15
		y16 = y-y15
		z16 = z-z15		
		
!  Q=eff*Volt*curr
		Q=0.85*29*270
!  a welding simulation 
		if(KSTEP.EQ.1) then

!     Goldak's volumetric heat source model for bead1 
				m=exp(-3*((x2)**2/(c)**2+(y2)**2/(bf)**2+(z2)**2/(c)**2))
				n=exp(-3*((x2)**2/(c)**2+(y2)**2/(br)**2+(z2)**2/(c)**2))
				FLUX(1)=(6*1.732*Q)*((1.4/bf)*m+(0.6/br)*n)/(pi*sqrt(pi)*(a*c))
!	for bead2	
				else
					if(KSTEP.EQ.3) then
				m=exp(-3*((x4)**2/(c)**2+(y4)**2/(bf)**2+(z4)**2/(c)**2))
				n=exp(-3*((x4)**2/(c)**2+(y4)**2/(br)**2+(z4)**2/(c)**2))
				FLUX(1)=(6*1.732*Q)*((1.4/bf)*m+(0.6/br)*n)/(pi*sqrt(pi)*(a*c))
!	for bead3				
					else
						if(KSTEP.EQ.5) then
					m=exp(-3*((x6)**2/(c)**2+(y6)**2/(bf)**2+(z6)**2/(c)**2))
					n=exp(-3*((x6)**2/(c)**2+(y6)**2/(br)**2+(z6)**2/(c)**2))
					FLUX(1)=(6*1.732*Q)*((1.4/bf)*m+(0.6/br)*n)/(pi*sqrt(pi)*(a*c))
!	for bead4				
						else
							if(KSTEP.EQ.7) then
						m=exp(-3*((x8)**2/(c)**2+(y8)**2/(bf)**2+(z8)**2/(c)**2))
						n=exp(-3*((x8)**2/(c)**2+(y8)**2/(br)**2+(z8)**2/(c)**2))
						FLUX(1)=(6*1.732*Q)*((1.4/bf)*m+(0.6/br)*n)/(pi*sqrt(pi)*(a*c))
				
!	for bead5				
							else
								if(KSTEP.EQ.9) then
							m=exp(-3*((y10)**2/(c)**2+(x10)**2/(bf)**2+(z10)**2/(c)**2))
							n=exp(-3*((y10)**2/(c)**2+(x10)**2/(br)**2+(z10)**2/(c)**2))
							FLUX(1)=(6*1.732*Q)*((1.4/bf)*m+(0.6/br)*n)/(pi*sqrt(pi)*(a*c))
				
!	for bead6				
								else
									if(KSTEP.EQ.11) then
								m=exp(-3*((y12)**2/(c)**2+(x12)**2/(bf)**2+(z12)**2/(c)**2))
								n=exp(-3*((y12)**2/(c)**2+(x12)**2/(br)**2+(z12)**2/(c)**2))
								FLUX(1)=(6*1.732*Q)*((1.4/bf)*m+(0.6/br)*n)/(pi*sqrt(pi)*(a*c))
!	for bead7			
									else
										if(KSTEP.EQ.13) then
									m=exp(-3*((y14)**2/(c)**2+(x14)**2/(bf)**2+(z14)**2/(c)**2))
									n=exp(-3*((y14)**2/(c)**2+(x14)**2/(br)**2+(z14)**2/(c)**2))
									FLUX(1)=(6*1.732*Q)*((1.4/bf)*m+(0.6/br)*n)/(pi*sqrt(pi)*(a*c))
!	for bead8				
										else
											if(KSTEP.EQ.15) then
										m=exp(-3*((y16)**2/(c)**2+(x16)**2/(bf)**2+(z16)**2/(c)**2))
										n=exp(-3*((y16)**2/(c)**2+(x16)**2/(br)**2+(z16)**2/(c)**2))
										FLUX(1)=(6*1.732*Q)*((1.4/bf)*m+(0.6/br)*n)/(pi*sqrt(pi)*(a*c))
										endif
									endif
								endif
							endif
						endif
					endif
				endif
			endif
		
	
	
	RETURN
	END