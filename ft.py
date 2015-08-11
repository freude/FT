import numpy as np
import matplotlib.pyplot as pl


def ft(t, f):

    dt = t[3] - t[2]

    # compute Fourier transform by numpy's FFT function
    g = np.fft.fft(f)
    # frequency normalization factor is 2*np.pi/dt
    w = np.fft.fftfreq(f.size)*2*np.pi/dt

    # In order to get a discretisation of the continuous Fourier transform
    # we need to multiply g by a phase factor
    g *= dt*np.exp(-1j*w*t[0])/(np.sqrt(2*np.pi))

    return (w, g)

def ft_corr(t, f):

    dt = t[3] - t[2]

    # frequency normalization factor is 2*np.pi/dt
    N = 256
    N = 1024
    w1 = np.fft.fftfreq(N)*2*np.pi/dt
    theta = w1*dt
    f1 = np.concatenate((f, np.zeros((N-f.size))))

    # compute Fourier transform by numpy's FFT function
    g = np.fft.fft(f1)
    cth=np.cos(theta)
    sth=np.sin(theta)
    ctth=cth**2-sth**2
    stth=2.0*sth*cth
    th2=theta*theta
    th4=th2*th2
    tmth2=3.0-th2
    spth2=6.0+th2
    sth4i=1.0/(6.0*th4)
    tth4i=2.0*sth4i
    corfac=tth4i*spth2*(3.0-4.0*cth+ctth)
    corfac[np.isnan(corfac)]=1.0
     
    a0r=sth4i*(-42.0+5.0*th2+spth2*(8.0*cth-ctth))
    a0i=sth4i*(theta*(-12.0+6.0*th2)+spth2*stth)
    a0r[np.isnan(a0r)]=-2.0/3.0
    a0i[np.isnan(a0i)]=0.0
    
    a1r=sth4i*(14.0*tmth2-7.0*spth2*cth)
    a1i=sth4i*(30.0*theta-5.0*spth2*sth)
    a1r[np.isnan(a1r)]=7.0/24.0
    a1i[np.isnan(a1i)]=0.0
    
    a2r=tth4i*(-4.0*tmth2+2.0*spth2*cth)
    a2i=tth4i*(-12.0*theta+2.0*spth2*sth)
    a2r[np.isnan(a2r)]=-1.0/6.0
    a2i[np.isnan(a2i)]=0.0
    
    a3r=sth4i*(2.0*tmth2-spth2*cth)
    a3i=sth4i*(6.0*theta-spth2*sth)
    a3r[np.isnan(a3r)]=1.0/24.0
    a3i[np.isnan(a3i)]=0.0

    cl=a0r*f[0]+a1r*f[1]+a2r*f[2]+a3r*f[3]
    sl=a0i*f[0]+a1i*f[1]+a2i*f[2]+a3i*f[3]
    cr=-a0r*f[-1]+a1r*f[-2]+a2r*f[-3]+a3r*f[-4]
    sr=-a0i*f[-1]-a1i*f[-2]-a2i*f[-3]-a3i*f[-4]

    arg=np.max(w1)*(t[-1]-t[0])
    c=np.cos(arg)
    s=np.sin(arg)
    corre=cl+c*cr-s*sr
    corim=sl+s*cr+c*sr
    cor = corre+1j*corim


    # In order to get a discretisation of the continuous Fourier transform
    # we need to multiply g by a phase factor
    g *= corfac
    g += cor
    print(corfac)
    g *= dt*np.exp(-1j*w1*t[0])/(np.sqrt(2*np.pi))

    return (w1, g)

if __name__ == "__main__":

    t=np.linspace(-110, 110, 450, endpoint=True)
    print(t)
    #Define function
    f=np.exp(-t**2)

    w,g = ft(t, f)
    w1,g1 = ft_corr(t, f)

    
    # Plot Result
    pl.scatter(w,np.real(g),color="r")
    pl.scatter(w1,np.real(g1),color="y")
    # For comparison we plot the analytical solution
    y = np.exp(-0.25*w**2)/np.sqrt(2.0)
    
    pl.plot(w,y,color="g")
    #pl.plot(w1,w1,color="g")

    pl.gca().set_xlim(-20, 20)
    pl.show()
    pl.close()
