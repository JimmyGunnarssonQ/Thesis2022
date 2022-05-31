import matplotlib.pyplot as plt 
from Functions import *


state = input("Declare type of graph:")


plt.grid(True)
if state == "mass":
    N1o, N2o, N4p = call(state)
    plt.plot(lambda3list, N4p, "--", color = "blue", label = r"$H_2^0$")
    plt.plot(lambda3list, N1o, color = "red", lw=4, label = r"$h$")
    plt.plot(lambda3list, N2o,  color = "blue", lw=4, label = r"$H$")
    plt.ylabel("Mass (GeV)")
    plt.title("Physical Higgs masses" +" as a function of " + r"$\lambda_3$")

elif state == "hVV":
    ang1p, ang2p = call(state)
    plt.plot(lambda3list, ang1p, lw=4, color = "red", label = r"$h$")
    plt.plot(lambda3list, ang2p, lw=4, color = "blue", label = r"$H$")
    plt.ylabel(r"$g_{hVV}^2$")
    plt.title("Vector boson couplings as a function of " + r"$\lambda_3$" + " with " + r"$t_\beta = \frac{1}{3}$")


else:
    print("False value")
plt.xlabel(r"$\lambda_3$")

plt.legend(title = "State:")
plt.show()
