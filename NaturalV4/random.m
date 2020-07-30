% Random Function used to create ux potential

%Creating and storing a random function

a0 = -2;
b0 = 2;
d = linspace(a0,b0,11);
well = rand(1,11);
sm = -5*rand/2;
sp = 5*rand/2;


save('random_potential.txt','d','well','-ascii')

