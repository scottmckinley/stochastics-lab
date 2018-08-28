




max_less_index_function = function(number, vector){
	return(max(which(vector <= number),0))
}


steppingOU_simulation = function(final_time, Delta, lambda, spring_const,sigma,simulation_time_step = 0.001, plot_it = FALSE){

	time_step = simulation_time_step;
	time_seq = seq(0, final_time, by = time_step)
	num_steps = length(time_seq)

	Tn = 0;	time_between_events = c()


	while(Tn < final_time){
		time_between_events = c(time_between_events, rexp(1, lambda))
		Tn = sum(time_between_events)		
	}
	time_events =cumsum(time_between_events)

	Zn = c(0,unlist(lapply(as.list(time_seq[-1]), max_less_index_function, vector = time_events )))

	
##=============================================
##Euler M
##=============================================
Xn_EM = c(rnorm(n = 1, mean = Zn[1],sd = 0.001) )
	
for(n in 1:(length(Zn)-1)){
	Xn_EM[n+1] = Xn_EM[n] - spring_const*(Xn_EM[n]- Zn[n])*time_step + sigma*sqrt(time_step)*rnorm(n = 1, mean =0, sd =1)
}
##=============================================
##Using the exact solution and the fact we know the motor position for 0 <= t <=T_final
##=============================================

Xn_exact = c(rnorm(n = 1, Zn[1], 0.001))
number_events = 0:length(time_between_events)

for( n in 1:(length(Zn)-1)){
	time_bds = c(time_seq[n], time_seq[n+1])
	multiple_event_index = intersect(which(time_events >= time_bds[1]), which(time_events <= time_bds[2]))
	
	if(length(multiple_event_index) == 0){ 
#Case 1: when no events happen between the interval [tn, t_n+1] so Z_n = Z_n+1
		dt_integral = Zn[n+1]*(1 - exp(-diff(time_bds)*spring_const))
	}else{
# Case 2: at least one event occurs between the interval [tn, t_n+1] so Z_n =! Z_n+1	
		time_event_n = sort(c(time_bds, time_events[multiple_event_index]))	
		interval_n  = diff(time_event_n)
		Z_continuous = c(number_events[multiple_event_index], Zn[n+1])
		dt_integral =sum(Z_continuous*(exp(-1*spring_const*(time_bds[2] - time_event_n)[-1]) - exp(-1*spring_const*(time_bds[2] - time_event_n)[-length(time_event_n)] )))
		}
	
	OU_var = (sigma^2/(2*spring_const))*(1 - exp(-spring_const*time_step*2))
	Xn_exact[n+1] = Xn_exact[n]*exp(-spring_const*time_step) + dt_integral + sqrt(OU_var)*rnorm(1, mean = 0, sd = 1)

	} #end of the simulate Xn - the cargo process


obs_index = seq(1, length(time_seq), by = Delta/time_step)
obs_time = time_seq[obs_index]
obs_Zn = Zn[obs_index]
obs_Xn_EM = Xn_EM[obs_index]; obs_Xn_exact= Xn_exact[obs_index]


if(plot_it == TRUE){
	
plot(time_seq, Zn, type = "l", xlab = "Time", ylab = "Position")
lines(time_seq, Xn_EM, col = "dodgerblue3")
lines(time_seq, Xn_exact, type = "l", col ="springgreen3")
	
plot(obs_time, obs_Zn, type = "l", xlab = "Time", ylab = "Position")
lines(obs_time, obs_Xn_EM, col = "dodgerblue3")
lines(obs_time, obs_Xn_exact, type = "l", col ="springgreen3")
legend("topleft", c("Motor", "EM", "Exact"), col = c("black", "dodgerblue3", "springgreen3"), lty = c(1,1,1))
}

return(list(time = obs_time, motor = obs_Zn, cargo = obs_Xn_exact, cargoEM = obs_Xn_EM))



}



# # final_time = 5; Delta = 0.1
# true_lambda = 0.81/0.008; spring_const = (1/100)*340/0.16; sigma = sqrt(2*0.0041/0.16)*500
# true_lambda = 0.81/0.008; spring_const = 100; sigma = 50

# my_sim  = steppingOU_simulation(final_time, Delta, true_lambda, spring_const, sigma, plot_it = TRUE)


