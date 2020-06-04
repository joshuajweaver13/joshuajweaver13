import matplotlib.pyplot as plt 
import math
from scipy.stats import norm
import numpy as np 
# Call price calculator

#  Establishing the input values
S0 = 150
σ = 0.30
r = 0.05
T = 5
K = 150

# Analyze the Change in Price when manipulating the Strike Price
Call_Prices_K = []
Put_Prices_K = []
for K in range(1, 400):
    # Calculate values d1 and d2
    d1 = float(((math.log((S0 / K)) + ((r + (0.5 * (σ ** 2))) * T)) / (σ * math.sqrt(T))))
    d2 = d1 - (σ * math.sqrt(T))

    # define a formula Φ(x) to find the Φ(.) values
    Φ = lambda x: norm.cdf(x)

    # Determine the call price 
    C0 = (S0 * Φ(d1)) - ((K * math.exp(-r * T)) * Φ(d2))
    Call_Prices_K.append(C0)

    # Deriving the Put price via Put-Call Parity
    P0 = C0 + (K * math.exp(-r * T)) - S0
    Put_Prices_K.append(P0)
plt.subplot(321)
plt.plot(Call_Prices_K, 'g-', label="Change in Call Price")
plt.plot(Put_Prices_K, 'b-', label='Change in Put Price')
plt.ylabel("Option Price")
plt.xlabel('Strike Price')
plt.title("Changing Strike Prices")

# Analyze the change in Price when manipulating the time from expiry
K = 150
Call_Prices_T = []
Put_Prices_T = []
for T in range(1,100):
    # Calculate values d1 and d2
    d1 = float(((math.log((S0 / K)) + ((r + (0.5 * (σ ** 2))) * T)) / (σ * math.sqrt(T))))
    d2 = d1 - (σ * math.sqrt(T))

    # define a formula Φ(x) to find the Φ(.) values
    Φ = lambda x: norm.cdf(x)

    # Determine the call price 
    C0 = (S0 * Φ(d1)) - ((K * math.exp(-r * T)) * Φ(d2))
    Call_Prices_T.append(C0)

    # Deriving the Put price via Put-Call Parity
    P0 = C0 + (K * math.exp(-r * T)) - S0
    Put_Prices_T.append(P0)
plt.subplot(322)
plt.plot(Call_Prices_T, 'g-')
plt.plot(Put_Prices_T, 'b-')
plt.ylabel("Option Price")
plt.xlabel('Time')
plt.title("Changing the time to expiry")

# Analyze the change in Price by maniplating the change in the c.c. risk free rate
T = 5
Call_Prices_r = []
Put_Prices_r = []
for r in np.arange(0.01, 0.150, 0.01):
    # Calculate values d1 and d2
    d1 = float(((math.log((S0 / K)) + ((r + (0.5 * (σ ** 2))) * T)) / (σ * math.sqrt(T))))
    d2 = d1 - (σ * math.sqrt(T))

    # define a formula Φ(x) to find the Φ(.) values
    Φ = lambda x: norm.cdf(x)

    # Determine the call price 
    C0 = (S0 * Φ(d1)) - ((K * math.exp(-r * T)) * Φ(d2))
    Call_Prices_r.append(C0)

    # Deriving the Put price via Put-Call Parity
    P0 = C0 + (K * math.exp(-r * T)) - S0
    Put_Prices_r.append(P0)
plt.subplot(323)
plt.plot(Call_Prices_r, 'g-')
plt.plot(Put_Prices_r, 'b-')
plt.ylabel("Option Price")
plt.xlabel('c.c. risk-free rate\n(in x%)')
plt.title("Changing c.c. risk-free rate")

r = 0.05
Call_Prices_σ = []
Put_Prices_σ = []
for σ in np.arange(0.01,1,0.01):
    # Calculate values d1 and d2
    d1 = float(((math.log((S0 / K)) + ((r + (0.5 * (σ ** 2))) * T)) / (σ * math.sqrt(T))))
    d2 = d1 - (σ * math.sqrt(T))

    # define a formula Φ(x) to find the Φ(.) values
    Φ = lambda x: norm.cdf(x)

    # Determine the call price 
    C0 = (S0 * Φ(d1)) - ((K * math.exp(-r * T)) * Φ(d2))
    Call_Prices_σ.append(C0)

    # Deriving the Put price via Put-Call Parity
    P0 = C0 + (K * math.exp(-r * T)) - S0
    Put_Prices_σ.append(P0)    
plt.subplot(324)
plt.plot(Call_Prices_σ, 'g-')
plt.plot(Put_Prices_σ, 'b-')
plt.ylabel("Option Price")
plt.xlabel('Volatility')
plt.title("Changing Volatility (σ)")

σ = 0.3
Call_Prices_S0 = []
Put_Prices_S0 = []
for S0 in range(1,251):
     # Calculate values d1 and d2
    d1 = float(((math.log((S0 / K)) + ((r + (0.5 * (σ ** 2))) * T)) / (σ * math.sqrt(T))))
    d2 = d1 - (σ * math.sqrt(T))

    # define a formula Φ(x) to find the Φ(.) values
    Φ = lambda x: norm.cdf(x)

    # Determine the call price 
    C0 = (S0 * Φ(d1)) - ((K * math.exp(-r * T)) * Φ(d2))
    Call_Prices_S0.append(C0)

    # Deriving the Put price via Put-Call Parity
    P0 = C0 + (K * math.exp(-r * T)) - S0
    Put_Prices_S0.append(P0)
plt.subplot(313)
plt.plot(Call_Prices_S0, 'g-', label="Change in Call Price")
plt.plot(Put_Prices_S0, 'b-', label='Change in Put Price')
plt.ylabel("Option Price")
plt.xlabel('Strike Price')
plt.title("Changing Underlying Stock's Spot Price")



plt.tight_layout(pad=1.0)
plt.subplots_adjust(top=0.875)
plt.subplots_adjust(bottom=0.175)
plt.suptitle("Manipulating the Variables within Black-Scholes",fontweight='heavy')
plt.figlegend(("European Call","European Put"),loc=8,ncol=2,fancybox=True,shadow=True)
plt.show()
