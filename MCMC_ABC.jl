########################################################################
## MCMC_ABC
## 05/01/2016
##
## 
########################################################################

function MCMC_ABC(initial, sd; iterations = 5, stepsize = 1, eps = 50, method = "AD", bucket_method = "NNB")
oecdf = observed_ecdf()
limit = maximum(train, 1)
p = length(initial)
chain = Matrix(iterations + 1, p + 1)
chain[1, p+1] = 0
chain[1, 1:p] = initial
rejection = 0
	for i = 1:iterations
		bar = Int64(round(i/iterations * 10))
		print("[ ", repeat("#", bar), repeat(" ", 10 - bar), " ]" ,"\u1b[1G")
		print(" ", (floor(i/iterations*10000))/100, "%", "\u1b[1G")
		print("[ ", repeat("#", bar), repeat(" ", 10 - bar), " ]", "\u1b[K")
		print(" ", (floor(i/iterations*10000))/100, "%", "\u1b[K")
		# print("iteration: ", i,"\n")
		proposal = proposalftn( chain[i, 1:p], sd )
		prob = exp( logPrior(proposal) - logPrior( chain[i,1:p] ) ) # this line is broken, term in expontential is 0 everytime
		u = rand( Uniform(0,1), 1 )[1]
		#@printf("rand is %2.4f \n", u)
		#@printf("prob is %2.4f \n", prob)
		if u < prob
			# print("theta proposal: ", proposal,"\n")
			u = rainfall_table_generator(proposal[4], proposal[5])
			rainfall_range = u[1] * 1000^2 * (3 * sqrt(3) / 2) * (2 / 3 * 0.1)^2 ;
			rainfall_table = u[2]
			yproposal = driver2(proposal[1:4], bucket_method, incidence, train, rainfall_range, heightdata)
			stat_j = zeros(length(yproposal)) 
			for j = 1:length(yproposal) # for each time period
				support = 0:stepsize:limit[j]
				if method == "AD" ## Anderson–Darling criterion
					weight = 1./(oecdf[j](support) ./ ( 1 - oecdf[j](support) ))
					stat_j[j] = sum( ( yproposal[j](support) - oecdf[j](support) ).^2 .* weight .* stepsize  )
				end
				if method == "CvM" ## Cramér–von Mises criterion
					stat_j[j] = sum( ( yproposal[j](support) - oecdf[j](support) ).^2 .*stepsize )
				end
				if method == "OdJ" ## Matching mean and variance
					weight = 1./(oecdf[j](support) ./ ( 1 - oecdf[j](support) ))
					mean_sim = sum( (1 - yproposal[j](support)).*stepsize )
					mean_obs = sum( (1 - oecdf[j](support)).*stepsize )
					var_sim = 2*sum( support.*yproposal[j](support).*stepsize ) - ( sum( (1 - yproposal[j](support)).*stepsize ) )^2
					var_obs = 2*sum( support.*oecdf[j](support).*stepsize ) - ( sum( (1 - oecdf[j](support)).*stepsize ) )^2
					# print("mean diff", 1e-3*(mean_sim - mean_obs)^2, "\n")
					# print("mean var", 1e-9*(var_sim - var_obs)^2, "\n")
					stat_j[j] = sum( ( yproposal[j](support) - oecdf[j](support) ).^2 .* weight .* stepsize  ) + 1e-3*(mean_sim - mean_obs)^2 + 1e-9*(var_sim - var_obs)^2
				end
			end
			stat = sum(stat_j)
			# print(method, " stat: ", stat, "\n")
			if stat < eps
				chain[i+1, 1:p] = proposal
				chain[i+1, 5] = stat
			else
				rejection += 1
				chain[i+1,:] = chain[i,:]
			end
		else
			rejection += 1
			chain[i+1,:] = chain[i,:]
		end
	end
ratio = (iterations - rejection)/iterations
avg = mapslices(mean,chain,1)
# print("\n")
# @printf("Ratio \t p \t tau \t scale \t a \t kappa \t delta \n")
# @printf("%2.2f \t %2.2f \t %2.2e \t %2.2f \t %2.2e \n",ratio, avg[1], avg[2],avg[3],avg[4],avg[5],avg[6] )
return chain
end