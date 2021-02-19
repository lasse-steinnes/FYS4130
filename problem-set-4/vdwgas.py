## programme to investigate properties of van der Waals gass
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

print("Press 1 to plot P vdw and P ideal")
print("Press 2 to plot T as func of V")
print("Press 3 to plot mu(P) and v(P)")
task = int(input("Choose task (int):"))

def p(tau, v):
    """
    reduced pressure
    as a function of reduced variables
    - tau: temp
    - v: volume
    """
    return (8*tau)/(3*v-1) - 3/(v**2)


def p_ideal(tau,v):
    """
    reduced pressure ideal gas
    as a function of reduced variables
    - tau: temp
    - v: volume
    """
    return 8*tau/(3*v)

def tau_vdw(p,v):
    """
    reduced temp
    """
    return (1/8)*(p*(3*v-1) + (9*v-3)/v**2)

def tau_ideal(p,v):
    """
    reduced temp ideal gas
    """
    return 3/8*v*p

def mu_ideal(tau,v):
    """
    mu for a particular tau as function of p
    """
    N = 10;
    mu = np.zeros(len(v)) + 0.1
    # calculate p
    p_id = p_ideal(tau,v)

    for i in range(len(p_id)-1):
        mu[i] += N*v[i]*(p_id[i+1]-p_id[i])
    return p_id, mu

def mu_vdw(tau,v):
    """
    mu for a particular tau
    """
    N = 10;
    mu = np.zeros(len(v)) + 0.1
    # calculate p
    p_vdw = p(tau,v)

    for i in range(len(p_vdw)-1):
        mu[i] += N*v[i]*(p_vdw[i+1]-p_vdw[i])
    return p_vdw, mu

if task == 1:
    ### make animation
    v = np.linspace(0.35,10,1000)

    # First set up the figure, the axis, and the plot element we want to animate
    fig,ax = plt.subplots(figsize= (6,6))
    fig.set_tight_layout(True)
    line, = ax.plot([],[], '-',label = "P reduced vdw")
    line2, = ax.plot([],[],'-', label = "P reduced ideal")
    ax.set_xlabel("$\hat{V}$")
    ax.set_ylabel("$\hat{P}$")
    ax.legend()

    # initialization function: plot the background of each frame
    def init():
        ax.set_xlim(0.0, 10)
        ax.set_ylim(-1, 10)
        return line,line2,

    # animation function.  This is called sequentially
    def animate(frame):
        #### using tau as frames
        ax.set_title("tau {:.2f}".format(frame))
        line.set_data(v, p(frame,v))
        line2.set_data(v,p_ideal(frame,v))
        return line, line2,

    anim = FuncAnimation(fig, animate, init_func=init,
                                   frames=np.linspace(0.0,1.5,100), interval=20)
    plt.show()

if task == 2:
        ### make animation
        v = np.linspace(0.1,8,1000)

        # First set up the figure, the axis, and the plot element we want to animate
        fig,ax = plt.subplots(figsize= (6,6))
        fig.set_tight_layout(True)
        line, = ax.plot([],[], '-',label = "T reduced vdw")
        line2, = ax.plot([],[],'-', label = "T reduced ideal")
        ax.set_xlabel("$\hat{V}$")
        ax.set_ylabel("$\hat{T}$")
        ax.legend()

        # initialization function: plot the background of each frame
        def init():
            ax.set_xlim(0.0, 8)
            ax.set_ylim(-1, 6)
            return line,line2,

        # animation function.  This is called sequentially
        def animate(frame):
            #### using p as frames
            ax.set_title("p {:.2f}".format(frame))
            line.set_data(v, tau_vdw(frame,v))
            line2.set_data(v,tau_ideal(frame,v))
            return line, line2,

        anim = FuncAnimation(fig, animate, init_func=init,
                                       frames=np.linspace(0.1,1.5,100), interval=20)
        plt.show()

if task == 3:
    ### make animation
    v = np.linspace(0.1,8,1000)

    # First set up the figure, the axis, and the plot element we want to animate
    fig,ax = plt.subplots(figsize= (6,6))
    fig.set_tight_layout(True)
    line, = ax.plot([],[], '-',label = "$\mu$ reduced vdw")
    line2, = ax.plot([],[],'-', label = "$\mu$ reduced ideal")
    ax.set_xlabel("$\hat{P}$")
    ax.set_ylabel("$\hat{\mu}$")
    ax.legend()

    # initialization function: plot the background of each frame
    def init():
        ax.set_xlim(0.0, 1.2)
        ax.set_ylim(-.2, .2)
        return line,line2,

    # animation function.  This is called sequentially
    def animate(frame):
        #### using tau as frames
        ax.set_title("tau {:.2f}".format(frame))
        pid,muid =  mu_ideal(frame,v)
        pvdw,muvdw = mu_vdw(frame,v)
        line.set_data(pvdw,muvdw)
        line2.set_data(pid, muid)
        return line, line2,

    anim = FuncAnimation(fig, animate, init_func=init,
                                   frames=np.linspace(0.8,1.2,50), interval=20)
    plt.show()

else:
     print("not an option")
