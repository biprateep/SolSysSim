"""

The program simulates our solar system and writes orbital parameters for a particular object in a file.
The program is written in Python 2.7 and needs a compatible installation of Visual Python to run.
It has been created by Sharba Bhattacharjee and Biprateep Dey at National Institute of Science Education and Research, Bhubaneswar, India.
The code is well commented but further explanation of the methods used can be found in the accompanying pdf document.
Any one is free to use the code in any way (s)he likes, but the authors take no responsibilty for any consequences due its usage.
Feel free to contact sharba.b@niser.ac.in or biprateep.d@niser.ac.in 

Last Modified: 12-12-2015
"""



import visual as v
import sys

G=6.67e-11 # Gravitational constant

# Defining the objects
sun=v.sphere(radius=3e10)
planet=[v.sphere(make_trail=True,retain=10000000) for i in range(10)] # Moon is also included here

# Data source: NASA Jet Propulsion Lab Horizons System
# Mass chart
sun.mass=1.989e30
planet[0].mass=3.302e23 # Mercury
planet[1].mass=4.8685e24 # Venus
planet[2].mass=5.97219e24 # Earth
planet[3].mass=6.4185e23 # Mars
planet[4].mass=1.89813e27 # Jupiter
planet[5].mass=5.68319e26 # Saturn
planet[6].mass=8.68103e25 # Uranus
planet[7].mass=1.0241e26 # Neptune
planet[8].mass=1.307e22 # Pluto
planet[9].mass=7.349e22 # Moon (Luna)

# Initial positions (1st Jan 2000, 12 a.m. California time)
planet[0].pos=v.vector(-2.105262111032039e10,-6.640663808353403e10,-3.492446023382954e9)
planet[1].pos=v.vector(-1.075055502695123e11,-3.366520720591562e9,6.159219802771119e9)
planet[2].pos=v.vector(-2.521092863852298e10,1.449279195712076e11,-6.164888475164771e5)
planet[3].pos=v.vector(2.079950549908331e11,-3.143009561106971e9,-5.178781160069674e9)
planet[4].pos=v.vector(5.989091594973032e11,4.391225931530510e11,-1.523254614945272e10)
planet[5].pos=v.vector(9.587063368200246e11,9.825652109121954e11,-5.522065682385063e10)
planet[6].pos=v.vector(2.158774703477132e12,-2.054825231595053e12,-3.562348723541665e10)
planet[7].pos=v.vector(2.514853420151505e12,-3.738847412364252e12,1.903947325211763e10)
planet[8].pos=v.vector(-1.477558339934327e12,-4.182460438550376e12,8.752694146056010e11)
planet[9].pos=v.vector(-2.552857888050620e10,1.446860363961675e11,3.593933517466486e7)

# Initial velocities
planet[0].vel=v.vector(3.66529870639384e4,-1.228983810111077e4,-4.368172898981951e3)
planet[1].vel=v.vector(8.891598046362434e2,-3.51592077412429e4,-5.318594054684045e2)
planet[2].vel=v.vector(-2.983983333368269e4,-5.207633918704476e3,6.169062303484907e-2)
planet[3].vel=v.vector(1.295003532851602e3,2.629442067068712e4,5.190097267545717e2)
planet[4].vel=v.vector(-7.901937610713569e3,1.116317695450082e4,1.306729070868444e2)
planet[5].vel=v.vector(-7.428885683466339e3,6.738814237717373e3,1.776643613880609e2)
planet[6].vel=v.vector(4.637648411798584e3,4.627192877193528e3,-4.285025663198061e1)
planet[7].vel=v.vector(4.465799984073191e3,3.075681163952201e3,-1.665654118310400e2)
planet[8].vel=v.vector(5.261925689692920e3,-2.648919644838698e3,-1.241873053678830e3)
planet[9].vel=v.vector(-2.927904627038706e4,-6.007566180814270e3,-1.577640655646029)

# The look
sun.color=(0.93,0.36,0.09)
sun.material=v.materials.emissive
planet[0].color=(0.45,0.45,0.45)
planet[1].color=(0.77,0.48,0.18)
planet[2].material=v.materials.earth
#planet[2].color=(0.18,0.18,0.4)
planet[3].color=(0.71,0.47,0.36)
planet[4].color=(0.65,0.53,0.42)
planet[5].color=(0.92,0.69,0.34)
planet[6].color=(0.73,0.87,0.90)
planet[7].color=(20,28,72)

sun.radius=35e9
planet[0].radius=3e9
planet[1].radius=6e9
planet[2].radius=6e9
planet[3].radius=3e9
planet[4].radius=20e9
planet[5].radius=15e9
planet[6].radius=10e9
planet[7].radius=10e9
planet[8].radius=3e9
planet[9].radius=1e9

# Planet for the moon
planet[9].mom=2 # Earth


t=0 # Initial time
dt=864*5 # Step size

for plnt in planet:
    plnt.newpos=v.vector(plnt.pos)

    plnt.doubpos=v.vector(plnt.pos) # For finding truncation error
    plnt.newdoubpos=v.vector(plnt.pos)
    plnt.doubvel=v.vector(plnt.vel)

def apair(i,j,odr): # Acceleration from interaction between planets
    return G*planet[j].mass*(planet[j].pos+odr*dt*planet[j].vel-planet[i].pos-odr*dt*planet[i].vel)/v.mag(planet[i].pos+odr*dt*planet[i].vel-planet[j].pos-odr*dt*planet[j].vel)**3

def av(i,odr): # Net acceleration
    a=-G*sun.mass*(planet[i].pos+odr*dt*planet[i].vel)/(v.mag(planet[i].pos+odr*dt*planet[i].vel))**3
    for j in range(len(planet)):
        if i!=j:
            a=a+apair(i,j,odr)
    return a

def doubapair(i,j,odr): # For truncation error
    return G*planet[j].mass*(planet[j].doubpos+odr*dt*planet[j].doubvel-planet[i].doubpos-odr*dt*planet[i].doubvel)/v.mag(planet[i].doubpos+odr*dt*planet[i].doubvel-planet[j].doubpos-odr*dt*planet[j].doubvel)**3

def doubav(i,odr):
    a=-G*sun.mass*(planet[i].doubpos+odr*dt*planet[i].doubvel)/(v.mag(planet[i].doubpos+odr*dt*planet[i].doubvel))**3
    for j in range(0,len(planet)):
        if i!=j:
            a=a+doubapair(i,j,odr)
    return a

# For planet data
ptf=int(raw_input("Which planet are you interested in (Mercury: 0,..., Pluto: 7) ")) # Planet to find (0: Mercury,..., 8 for Pluto)
times=float(raw_input("How many data points do you want for the planet: ")) # Number of times data are taken

timc=0 # Counter

year=[] # List of orbital periods
per=[] # List of perihelion distances
aph=[] # List of aphelion distances

pos0=v.vector(planet[ptf].pos) # Initial position of chosen planet

pe1=v.mag(pos0) # Variable for perihilion distance
ap1=v.mag(pos0) # Variable for aphilion distance
dirc0=v.dot(pos0,pos0) # Dot product (changes sign when angle crosses pi/2: keeping track of revolutions)

t_init=0 # Initial time for starting calculation (when the planet is perpendicular to initial position vector for the first time)


# For moon data
timoon=float(raw_input("How many data points do you want for the moon: ")) # Number of times moon data are taken
#mtf=int(raw_input("Which moon are you interested in (Moon: 0, Deimos: 1) ")) # Moon to find
mtf=9 # Moon to find: only one anyway
timcm=0 # Counter
month=[] # List of lunar orbital periods
posm0=planet[mtf].pos-planet[planet[mtf].mom].pos # Initial position of moon relative to its planet
dircm0=v.dot(posm0,posm0) # Dot product
tm_init=0 # Initial time for taking moon data
peg=[] # List of prigee distances
apg=[] # List of apogee distance
pg1=v.mag(posm0) # Variable for perigee distance
ag1=v.mag(posm0) # Variable for apogee distance

iteration_count=0

Roundp=[]
Truncp=[]
Roundm=[]
Truncm=[]

# Loop for orbital motion
while True: # Infinite loop
    v.rate(20000000000000000000000)
    for i in range(len(planet)): # Second order Runge-Kutta method
        planet[i].newpos=planet[i].pos+planet[i].vel*dt+0.5*av(i,0)*dt*dt
        planet[i].vel=planet[i].vel+dt*av(i,0.5)
    for plnt in planet:
        plnt.pos=v.vector(plnt.newpos)
    t=t+dt
    planet[2].rotate(angle=2*dt*v.pi/86400,axis=v.vector(0,0,1)) # Earth's rotation
    
    if len(year)<times or len(month)<timoon: # For truncation error
        if iteration_count%2 != 0:
            for i in range(len(planet)):
                planet[i].newdoubpos=planet[i].doubpos+planet[i].doubvel*2*dt+0.5*av(i,0)*dt*dt*4 # For truncation error
                planet[i].doubvel=planet[i].doubvel+2*dt*av(i,1)
            for plnt in planet:
                plnt.doubpos=v.vector(plnt.newdoubpos)
    
    # To find data for planet        
    if len(year)<times: # Data to be taken only for given number of times
        dirc=v.dot(planet[ptf].pos,pos0) # Dot product between initial and current position vectors of planet
        pe1=min(pe1,v.mag(planet[ptf].pos)) # Will give perihelion distance over each revolution
        ap1=max(ap1,v.mag(planet[ptf].pos)) # Will give aphelion distance over each revolution
        if dirc*dirc0<=0: # Sign changes for perpendicular direction; twice in one revolution
            timc+=1
            if timc==1:
                t_init=t # First time it becomes perpendicular; orbital period measured starting from this time
            if timc%2 != 0 and timc>1: # One revolution after t_init
                year.append((t-t_init)/86400.0-sum(year)) # Orbital period
                per.append(pe1/1000) # Perihelion distance
                aph.append(ap1/1000) # Aphelion distance
                Truncp.append(v.mag(planet[ptf].pos-planet[ptf].doubpos)/7000) # Truncation error
                Roundp.append((9000*iteration_count)**0.5*sys.float_info.epsilon) # Round-off error
                pe1=v.mag(planet[ptf].pos) # Initializing for next perihelion
                ap1=v.mag(planet[ptf].pos) # Initializing for next aphelion
                if (len(year)==times): # Required number of data points have been taken for planet
                    f=open('Planet_data.txt','w') # Output data to file
                    f.write('DATA OF PLANET ' + str(ptf) + '\nSl. no\tTime period\tPerihelion distance\tAphelion distance\tTruncation error\tRound-off error\n\t(earth day)\t(km)\t\t\t(km)\n')
                    for i in range(len(year)):
                        f.write(repr(i+1)+'\t'+str(year[i])+'\t\t'+str(per[i])+'\t\t'+str(aph[i])+'\t\t'+str(Truncp[i])+'\t\t'+str(Roundp[i])+'\n')
                    f.write('Average\t'+str(sum(year)/len(year))+'\t\t'+str(sum(per)/len(per))+'\t\t'+str(sum(aph)/len(aph)))
                    f.close()
                    print "Time period of planet (earth day)=",sum(year)/len(year) # Output average data on screen
                    print "Perihelion distance (km)=",'%e' % (sum(per)/len(per))
                    print "Aphelion distance (km)=",'%e' %  (sum(aph)/len(aph))
        dirc0=dirc
        
    # To find data for moon (Exactly similar logic)
    if len(month)<timoon: # Take data only for given number of times
        dircm=v.dot(planet[mtf].pos-planet[planet[mtf].mom].pos,pos0) # Dot product
        pg1=min(pg1,v.mag(planet[mtf].pos-planet[planet[mtf].mom].pos)) # For perigee
        ag1=max(ag1,v.mag(planet[mtf].pos-planet[planet[mtf].mom].pos)) # For apogee
        if dircm*dircm0<=0: # Crosses pi/2
            timcm+=1 # Half revolution
            if timcm==1:
                tm_init=t # Initial time
            if timcm%2!=0 and timcm>1:
                month.append((t-tm_init)/86400.0-sum(month)) # Orbital period of moon
                peg.append(pg1/1000) # Perigee
                apg.append(ag1/1000) # Apogee
                Truncm.append(v.mag(planet[mtf].pos-planet[mtf].doubpos)/7000) # Truncation error
                Roundm.append((9000*iteration_count)**0.5*sys.float_info.epsilon) # Round-off error
                pg1=v.mag(planet[mtf].pos-planet[planet[mtf].mom].pos) # Initializing for next perigee
                ag1=v.mag(planet[mtf].pos-planet[planet[mtf].mom].pos) # Initializing for next apogee
                if (len(month)==timoon): # Required number of data points have been taken for moon
                    fm=open('Moon_data.txt','w') # Output data to file
                    fm.write('DATA OF MOON\nSl. no\tTime period\tPerigee distance\tApogee distance\tTruncation error\tRound-off error\n\t(earth day)\t(km)\t\t\t(km)\n')
                    for i in range(len(month)):
                        fm.write(repr(i+1)+'\t'+str(month[i])+'\t\t'+str(peg[i])+'\t\t'+str(apg[i])+'\t'+str(Truncm[i])+'\t\t'+str(Roundm[i])+'\n')
                    fm.write('Average\t'+str(sum(month)/len(month))+'\t\t'+str(sum(peg)/len(peg))+'\t\t'+str(sum(apg)/len(apg)))
                    fm.close()
                    print "Time period for moon (earth day)=",sum(month)/len(month) # Output average data on screen
                    print "Perigee distance (km)=",'%e' % (sum(peg)/len(peg))
                    print "Apogee distance (km)=",'%e' %  (sum(apg)/len(apg))
        dircm0=dircm

    iteration_count+=1
