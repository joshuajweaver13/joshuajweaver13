import matplotlib.pyplot as plt 
import math
from scipy.stats import norm

# Call price calculator

#  Establishing the input values
S0 = 160
σ = 0.30
r = 0.1
K = 150
T = 2

# Calculate values d1 and d2
d1 = float(((math.log((S0 / K)) + ((r + (0.5 * (σ ** 2))) * T)) / (σ * math.sqrt(T))))
d2 = d1 - (σ * math.sqrt(T))

# define a formula Φ(x) to find the Φ(.) values
Φ = lambda x: norm.cdf(x)

# Determine the call price 
C0 = (S0 * Φ(d1)) - ((K * math.exp(-r * T)) * Φ(d2))
C0_formatted = f'{(C0):,.2f}'

# Deriving the Put price via Put-Call Parity
P0 = C0 + (K * math.exp(-r * T)) - S0
P0_formatted = f'{(P0):,.2f}'

# Creating the Call Payoff and profits to be plot in matplotlib
C0_Payoff = []
C0_Profit = []
for St in range(0, S0 * 2):
    x = max(St - K, 0)
    C0_Payoff.append(x)
    y = x - C0
    C0_Profit.append(y)
plt.subplot(121)
plt.plot(C0_Payoff,'b--', label="Payoff to the European Call Option")
plt.plot(C0_Profit, 'b-', label="Profit from the European Call Option")
plt.title("European Call Option")
plt.ylabel('Payoff')
plt.xlabel('Spot Price')

#  Creating the Put Payoff and profits to be plotted in matplotlib
P0_Payoff = []
P0_Profit = []
for St in range(0, S0 * 2):
    x = max(K - St, 0)
    P0_Payoff.append(x)
    y = x - P0
    P0_Profit.append(y)
plt.subplot(122)
plt.plot(P0_Payoff,'g--', label="Payoff to the European Put Option")
plt.plot(P0_Profit, 'g-', label="Profit from the European Put Option")
plt.title("European Put Option")
plt.ylabel('Payoff')
plt.xlabel('Spot Price')


# Formatting the programs outputs
print("The Price of the European Call option is $" + C0_formatted)
print()
print("Assuming No-Arbitrage Pricing and Put-Call Parity:")
print("The Price of the European Put option is $" + P0_formatted)


plt.tight_layout(pad=5.0)
plt.figlegend(loc=8,ncol=2,fancybox=True,shadow=True,framealpha=0.9,borderaxespad=0.9)
plt.suptitle("European Option Profit and Payoffs\nBased on the Black-Scholes Pricing Model")
plt.show()
