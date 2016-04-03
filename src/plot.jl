#####################################
#
#	Plotting functions
#
#####################################

function plot_ecdf(period, heightdata, runoff_ranking, rainfall_range; p = 1, tau = 100, a = 0.3, b = 6)
   infiltration = infiltration_gen(heightdata, 1/0.0056/20*1000^2*2/3/sqrt(3)/b);
   runoff_directions = direction_simulation2(heightdata, p);
   lookup_table = lookup_table_sim(tau, a, runoff_directions, infiltration, rainfall_range)
   final_fantasy = lookup_table * rainfall_table * 1e-6
   final_fantasy_II = missingdata(final_fantasy)
   emp_cdf = ecdf(final_fantasy_II[:,period])
   obs_cdf = observed_ecdf()[period]
   lim = maximum(train[:,period]) + 1
   figure()
   plot(emp_cdf(0:0.1:lim), color = "red", label = "simulation")
   plot(obs_cdf(0:0.1:lim), color = "black", label = "observed")
   legend(loc = 4)
   xlabel("Runoff")
end

function plot_all_ecdf(heightdata, runoff_ranking, rainfall_range; p = 1, tau = 100, a = 0.3, b = 6)
   infiltration = infiltration_gen(heightdata, 1/0.0056/20*1000^2*2/3/sqrt(3)/b);
   runoff_directions = direction_simulation2(heightdata, p);
   lookup_table = lookup_table_sim(tau, a, runoff_directions, infiltration, rainfall_range)
   final_fantasy = lookup_table * rainfall_table * 1e-6
   final_fantasy_II = missingdata(final_fantasy)
   lim = maximum(train)
   figure()
   for period = 1:23
      emp_cdf = ecdf(final_fantasy_II[:,period])
      obs_cdf = observed_ecdf()[period]
      plot(emp_cdf(0:0.1:lim), color = "red", label = "simulation", alpha = 0.5)
      plot(obs_cdf(0:0.1:lim), color = "black", label = "observed", alpha = 0.5)
   end
   xlabel("Runoff")
end

function plot_hist(period, heightdata, runoff_ranking, rainfall_range; p = 1, tau = 100, a = 0.3, b = 6)
   infiltration = infiltration_gen(heightdata, 1/0.0056/20*1000^2*2/3/sqrt(3)/b);
   runoff_directions = direction_simulation2(heightdata, p);
   lookup_table = lookup_table_sim(tau, a, runoff_directions, infiltration, rainfall_range)
   final_fantasy = lookup_table * rainfall_table * 1e-6
   final_fantasy_II = missingdata(final_fantasy)
   figure()
   plt[:hist](final_fantasy_II[:,period], bins = minimum(final_fantasy_II[:,period]):20:(maximum(final_fantasy_II[:,period])+50), alpha = 0.5, label = "simulated")
   plt[:hist](train[:,period], bins = minimum(final_fantasy_II[:,period]):20:(maximum(final_fantasy_II[:,period]) + 50), alpha = 0.5, label = "observed")
   xlabel("Runoff")
   ylabel("Frequency")
   legend(loc = 1)
end

function plot_all_hist(heightdata, runoff_ranking, rainfall_range; p = 1, tau = 100, a = 0.3, b = 6)
   infiltration = infiltration_gen(heightdata, 1/0.0056/20*1000^2*2/3/sqrt(3)/b);
   runoff_directions = direction_simulation2(heightdata, p);
   lookup_table = lookup_table_sim(tau, a, runoff_directions, infiltration, rainfall_range)
   final_fantasy = lookup_table * rainfall_table * 1e-6
   final_fantasy_II = missingdata(final_fantasy)
   figure(23)
   subplot(211)
   plt[:hist](final_fantasy_II)
   title("Simulation")
   xlabel("Runoff")
   ylabel("Frequency")
   subplot(212)
   plt[:hist](train)
   title("Observed")
   xlabel("Runoff")
   ylabel("Frequency")
end

function plot_correlogram(chain, lag = 20)
   figure()
   subplot(221)
      plt[:acorr](chain[:, 1], maxlags = lag)
      ylabel("p")
   subplot(222)
      plt[:acorr](chain[:, 2], maxlags = lag)
      ylabel("tau")
   subplot(223)
      plt[:acorr](chain[:, 3], maxlags = lag)
      ylabel("scale")
   subplot(224)
      plt[:acorr](chain[:, 4], maxlags = lag)
      ylabel("a")
end

function plot_data(theta, method = "WAVG")
	infiltration = infiltration_gen(heightdata, theta[4])
	runoff_directions = direction_simulation3(heightdata, incidence, theta[1]);
	lookup_table = lookup_table_sim2(theta[2], theta[3], method, incidence, runoff_directions, ranked, infiltration, rainfall_range)
	final_fantasy = lookup_table * rainfall_table * 1e-6 ######## mm3 to L
	final_fantasy_II = missingdata(final_fantasy)
	diff = train - final_fantasy_II
	figure()
	subplot(2,2,1)
		pcolormesh(diff)
		colorbar()
		ylabel("diff")
	subplot(2,2,2)
	pcolormesh(train)
		colorbar()
		ylabel("train")
	subplot(2,2,3)
	pcolormesh(final_fantasy_II)
		colorbar()
		ylabel("sim")
end

function plot_LUT(heightdata, runoff_ranking, rainfall_range; p = 1, tau = 100, a = 0.3, b = 6, color = "RdPu")
infiltration = infiltration_gen(heightdata, 1/0.0056/20*1000^2*2/3/sqrt(3)/b);
runoff_directions = direction_simulation2(heightdata, p);
lookup_table = lookup_table_sim(tau, a, runoff_directions, infiltration, rainfall_range)
figure()
pcolormesh(lookup_table, cmap = ColorMap(color))
colorbar()
axis([0, size(lookup_table,2), 0, size(lookup_table,1)])
xlabel("Rainfall level")
ylabel("Bucket #")
end

function plot_output(heightdata, runoff_ranking, rainfall_range; p = 1, tau = 100, a = 0.3, b = 6, color = "RdPu")   
   infiltration = infiltration_gen(heightdata, 1/0.0056/20*1000^2*2/3/sqrt(3)/b);
   runoff_directions = direction_simulation2(heightdata, p);
   lookup_table = lookup_table_sim(tau, a, runoff_directions, infiltration, rainfall_range)
   final_fantasy = lookup_table * rainfall_table * 1e-6
   final_fantasy_II = missingdata(final_fantasy)
   figure()
   pcolormesh(final_fantasy_II, cmap = color)
   colorbar()
   axis([0, size(lookup_table,2), 0, size(lookup_table,1)])
   xlabel("period")
   ylabel("bucket #")
end

function plot_posterior(chain)
   figure()
   subplot(221)
      plt[:hist](chain[:,1], alpha = 0.6, bins = 30, label = "p posterior")
      xlabel("p")
   subplot(222)
      plt[:hist](chain[:,2], alpha = 0.6, bins = 30, label = "tau posterior")
      xlabel("tau")
   subplot(223)
      plt[:hist](chain[:,3], alpha = 0.6, bins = 30, label = "scale posterior")
      xlabel("scale")
   subplot(224)
      plt[:hist](chain[:,4], alpha = 0.6, bins = 30, label = "a posterior")
      xlabel("a")
end

function plot_runoff(heightdata, runoff_ranking, rainfall_range; plot_runoff_one = 1, contour_if_one = 1, rainfall = 80000, p = 1, tau = 100, a = 0.3, b = 6, color = "RdPu")
infiltration = infiltration_gen(heightdata, b);
runoff_directions = direction_simulation2(heightdata, p);
runoff_sim_output = runoff_sim2(runoff_ranking, rainfall, infiltration, incidence, runoff_directions, tau, a);
x = size(runoff_sim_output,2)
y = size(runoff_sim_output,1)
if plot_runoff_one == 1
	figure()
	p = zeroreplacer(runoff_sim_output)
	pcolormesh(log(1.5,1+p), cmap = ColorMap(color))#, vmin = 20, vmax = 70 )
	colorbar()
	axis([0, x , 0, y])
	title("log transformed runoff")
	plt[:scatter](bucket[:,2], bucket[:,1]) 
	if contour_if_one == 1
		CS = plt[:contour](heightdata2, colors = "white", alpha = 1)
		clabel(CS, inline = 1, fontsize = 10)
	end
end
if plot_runoff_one == 2
	figure()
	p = zeroreplacer(runoff_sim_output)
	pcolormesh(p, cmap = ColorMap(color))#, vmin = 20, vmax = 70 )
	colorbar()
	axis([0, x , 0, y])
	title("runoff")
	plt[:scatter](bucket[:,2], bucket[:,1]) 
	if contour_if_one == 1
		CS = plt[:contour](heightdata2, colors = "white", alpha = 1)
		clabel(CS, inline = 1, fontsize = 10)
	end
end
if plot_runoff_one == 3
	for i in 1:23
	runoff_sim_output = runoff_sim2(runoff_ranking, rainfall_range[i], infiltration, incidence, runoff_directions, tau, a);
	p =  zeroreplacer(runoff_sim_output)
	pcolormesh(log(1.5,1+p), cmap = ColorMap(color), vmin = 20, vmax = 70 ) 
	colorbar()
	axis([0, x , 0, y])
	if contour_if_one == 1
		CS = plt[:contour](heightdata2, colors = "white", alpha = 1)
		clabel(CS, inline = 1, fontsize = 10)
	end
	savefig(string("rainfall level ", i,".png"), dpi = 72)
	close()
	end
end
end


function plot_runoff_movie(logscale = 0; save = 0, interval = 100, dpi = 250, bitrate = 1800, p = 1, tau = 100, a = 0.3, b = 6 ,color = "viridis")
infiltration = infiltration_gen(heightdata, 1/0.0056/20*1000^2*2/3/sqrt(3)/b);
runoff_directions = direction_simulation2(heightdata, p);
x = size(runoff_directions,2)
y = size(runoff_directions,1)
	ims = []
	fig = figure()
	for i in 1:23
	runoff_sim_output = runoff_sim2(runoff_ranking, rainfall_range[i], infiltration, incidence, runoff_directions, tau, a);
	p =  zeroreplacer(runoff_sim_output)
	if logscale == 0
		im = imshow( p, cmap = ColorMap(color), 
			vmin = 0, 
			vmax = 1.5e10, 
			origin = "lower", 
			aspect = "auto" ) 
	elseif logscale == 1
		im = imshow( log(1.5, 1 + p) , cmap = ColorMap(color), 
			vmin = 0, 
			vmax = 70, 
			origin = "lower", 
			aspect = "auto" ) 
	end
	if i == 1
		colorbar()
		if logscale == 0
			CS = plt[:contour](heightdata2, colors = "white", alpha = 0.1)
		elseif logscale == 1
			CS = plt[:contour](heightdata2, colors = "white", alpha = 0.3)
		end
		clabel(CS, inline = 1, fontsize = 10)
	end
	push!(ims, PyCall.PyObject[im])
	end	
	anim = animation.ArtistAnimation(fig, ims, interval = interval, blit = false, repeat = true)
	if save == 1
		anim[:save]("RAINANIMATION.mp4", dpi = dpi, bitrate = bitrate) 
	end
end

function plot_trace(chain)
   figure()
   subplot(221)
      plot(chain[:,1], label = "p posterior")
      ylabel("p")
   subplot(222)
      plot(chain[:,2], label = "tau posterior")
      ylabel("tau")
   subplot(223)
      plot(chain[:,3], label = "scale posterior")
      ylabel("scale")
   subplot(224)
      plot(chain[:,4], label = "a posterior")
      ylabel("a")
end


function qqplot(theta, period, logscale = 1)
if period > 23
	print("Period must be limited from 1 to 23")
else
	infiltration = infiltration_gen(heightdata, 173/20 * 1000^2 * (3 * sqrt(3) / 2) * (2 / 3 * 0.1)^2 /theta[4])
	runoff_directions = direction_simulation2(heightdata, theta[1]);
	lookup_table = lookup_table_sim(theta[2], theta[3], runoff_directions, infiltration, rainfall_range)
	final_fantasy = lookup_table * rainfall_table * 1e-6 ######## mm3 to L
	final_fantasy_II = missingdata(final_fantasy)
	x_ord = sort(final_fantasy_II[:,period])
	y_ord = sort(train[:,period])
	if logscale == 0
		figure()
		plot(x_ord, y_ord, 
			linewidth = 0, 
			marker = "o",  
			markersize = 7,
			markeredgewidth = 1,
			markeredgecolor = "black",
			markerfacecolor = "none")
		xlabel("Simulation")
		ylabel("Observed")
		title("Normal scale")
		x = linspace(0, max(x_ord[34], y_ord[34]), 50)
		plot(x, x, ls = "dashed", color = "black")
	elseif logscale == 1
		figure()
		plot(log(x_ord), log(y_ord), 
			linewidth = 0, 
			marker = "o",  
			markersize = 7,
			markeredgewidth = 1,
			markeredgecolor = "black",
			markerfacecolor = "none")
		xlabel("log(Simulation)")
		ylabel("log(Observed)")
		title("Log scale")
		x = linspace(0, log(1 + max(x_ord[34], y_ord[34])), 50)
		plot(x, x, ls = "dashed", color = "black")
	end
end
end

function qqplot_all(theta, method)
	infiltration = infiltration_gen(heightdata, theta[4])
	@time runoff_directions = direction_simulation2(heightdata, incidence, theta[1]);
	@time lookup_table = lookup_table_sim(theta[2], theta[3], method, incidence, runoff_directions, ranked, infiltration, rainfall_range);
		final_fantasy = lookup_table * rainfall_table * 1e-6 ######## mm3 to L
		final_fantasy_II = missingdata(final_fantasy)
		x_ord = sort(final_fantasy_II,1)
		y_ord = sort(train,1)
		figure()
		for i = 1:23 # 23 periods
			subplot(4,6,i)
				plot(log(50000,1 + x_ord[:,i]), log(50000,1 + y_ord[:,i]), 
				linewidth = 0, 
				marker = "o",  
				markersize = 5,
				markeredgewidth = 0.5,
				markeredgecolor = "black",
				markerfacecolor = "none")
				x = linspace(0, log(50000,1 + max(x_ord[34,i], y_ord[34,i])), 50)
				plot(x, x, ls = "dashed", color = "black")
				title(string("period ", i))
				xlabel("sim")
				ylabel("obs")
		end
end

