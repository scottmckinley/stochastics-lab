final_T = 10
num_steps = 1000
D = 1

t = seq(0,final_T,length = num_steps)
dt = t[2]-t[1]
noise = rnorm(num_steps,0,1)

x = sqrt(dt*D)*cumsum(noise)

plot(t,x,type="s",ylim = c(-2*sqrt(final_T*D),2*sqrt(final_T*D)),
     main="Brownian Motion")