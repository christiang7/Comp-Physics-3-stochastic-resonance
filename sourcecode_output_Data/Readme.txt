execution:

g++ Main.cpp -o stoch_reso && ./stoch_reso 20 1 1 0.05 10000 0

description of the parameters:

./stoch_reso p1 p2 p3 p4 p5 p6

p1 - T - duration of simulation
p2 - A - A, amplitude of the oscillator
p3 - omega - omega, oscillation of the oscillatior
p4 - sigma - sigma, amplitude of the noise
p5 - N - number of calculations
p6 - x_0, start value of x

history:

29.1.2019:
    new parameters:
        N - number of calculations 
        T - duration of simulation
    delta_t is calculated with the formula delta_t = T/N
    the white noise calculated with mu = 0 and variance = sqrt(tstep)
