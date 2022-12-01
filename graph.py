import matplotlib.pyplot as plt
import numpy as np

def from_file(str):
    f, l = open(str, "r"), []
    for line in f:   
        for w in line.split():         
            l.append(float(w))
    f.close()
    return l

def main():
    l0 = from_file("./res/0.txt")
    l1 = from_file("./res/1.txt")
    l2 = from_file("./res/2.txt")
    l3 = from_file("./res/3.txt")
    l4 = from_file("./res/4.txt")
    l5 = from_file("./res/5.txt")

    fig, ax = plt.subplots()

    ax.plot(l0, label='Initial condition')
    ax.plot(l1, label='Godunov scheme')
    ax.plot(l2, label='MacCormack scheme')
    ax.plot(l3, label='Holodnov scheme')
    ax.plot(l4, label='Cubic approximation scheme')
    ax.plot(l5, label='Hybrid scheme')
    ax.legend()

    fig.set_figheight(5)
    fig.set_figwidth(8)
    plt.show()
    
if __name__ == "__main__":
    main()