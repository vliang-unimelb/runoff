

function filler(x)
	for i = 1:size(x,1)
		if i % 2 == 1
			x[i,:] = repmat([1,0], div( size(x,2), 2 ) )
		else
			x[i,:] = repmat([0,1], div( size(x,2), 2 ) )
		end
	end
	return x
end

function hex_extract(data)
  x = zeros(size(data,1), size(data,2))
  if size(data,2) % 2 == 1
    println("ODD COLUMNS")
    x = hcat(x, zeros(size(data,1)))
    xfill = filler(x)
    x = x[ : ,1:( size(xfill,2) -1 ) ] 
    return x
  else 
    println("EVEN COLUMNS")
    return filler(x)
  end
end


function findneibheight(pos, all_height)
neib = zeros(6)
    neib[1] = all_height[ pos[1] - 2, pos[2] ] # N
    neib[2] = all_height[ pos[1] - 1, pos[2] + 1 ]  # NE
    neib[3] = all_height[ pos[1] + 1, pos[2] + 1 ]  # SE
    neib[4] = all_height[ pos[1] + 2, pos[2] ]  # S
    neib[5] = all_height[ pos[1] + 1, pos[2] - 1 ] # SW
    neib[6] = all_height[ pos[1] - 1, pos[2] - 1 ]  # NW
  return neib
end

function ranking(data)
heightdata = copy(data)
n1 = size(heightdata,1)
n2 = size(heightdata,2)
heightdata[1:2,:] = 0 # top border
heightdata[:,1:2] = 0 #left
heightdata[:, (n2 - 1) : n2] = 0 # right
heightdata[(n1 - 1) : n1, : ] = 0 # bottom
ranked = Array(Int64, 0)
sorted = sortperm(vec(heightdata), rev=true)[ 1: div((n1 - 4) * (n2 - 4), 2) ]
for j in 1:length(sorted)
	index = ind2sub( size(heightdata), sorted[j] )
	ranked = push!(ranked, index[1], index[2] )
end
ranked = transpose( reshape( ranked, ( 2, div( length(ranked), 2) ) ) )
return ranked
end


function infiltration_gen(data, scale)
rho = 173/20 * 1000^2 * (3 * sqrt(3) / 2) * (2 / 3 * 0.1)^2 / scale
e = Exponential(rho)
evalue = rand(e,length(data))
infiltration = reshape(evalue, (size(data,1), size(data,2)) )
return infiltration
end


function killnegative(die)
	for i = 1:6
		if die[i] < 0
			die[i] = 0
		end
	end
	return die
end	


function findpositive(live)
	r = Array(Int64, 0)
	for i = 1:6
		if live[i] > 0
			push!(r, i)
		end
	end
	return r
end


function my_all(diff) 
	for i = 1:6
		if diff[i] > 0
		return false
		end
	end
	return true
end 


function prob_gen(diff,p)
	r = Vector(0)
	@simd for i = 1:length(diff)
		if diff[i] > 0
			push!(r,diff[i]^p)
		end
	end
	return map(Float64,r/sum(r))
end



function direction_simulation2(heightdata, incidence, p = 1)
n1 = size(heightdata,1)
n2 = size(heightdata,2)
check = copy(incidence)
check[1:2,:] = 0 # top border
check[:,1:2] = 0 #left
check[:, (n2 - 1) : n2] = 0 # right
check[(n1 - 1) : n1, : ] = 0 # bottom
check = find(check .== 1)
d1_runoff = zeros(n1,n2)
d2_runoff = zeros(n1,n2)
	for z = 1:length(check)
		coord = ind2sub(incidence, check[z])
		i = coord[1]
		j = coord[2]
		neib_height = findneibheight([i, j], heightdata)
		diff = heightdata[i, j] - neib_height
		diff = killnegative(diff)
		if my_all(diff)
		else
			m = zeros(6)
			want = findpositive(diff)
			prob = prob_gen(diff, p)
			index = findfirst(rand( Multinomial( 1, prob ), 1 ), 1)
			d1 = want[index] # alteratively can do matrix multiplication
			d1_runoff[i,j] = d1
			if length(want) > 1
				deleteat!(want, index)
				deleteat!(diff, d1)
				prob = prob_gen(diff, p)
				d2 = want[findfirst(rand( Multinomial( 1, prob ), 1 ), 1)]
			else
				d2 = d1
			end		
			d2_runoff[i,j] = d2
		end	
	end
return typeof(d1_runoff)[d1_runoff, d2_runoff]
end


function runoff_sim2(ranked, rainfall, infiltration, incidence, runoff_directions, tau, a=1/4)
d1_direction = runoff_directions[1]
d2_direction = runoff_directions[2]
runoff = (rainfall - infiltration).*incidence
siz = size(ranked, 1)
	for k = 1 : siz 
		i = ranked[k, 1] 
		j = ranked[k, 2]
		d1 = 0
		d2 = 0
		runoff[i,j] = max( runoff[i,j], 0 )
		runoffij = runoff[i,j]
		if runoffij < tau
			d1 = runoffij
			d2 = 0
		else
			d1 = tau + (runoffij - tau) * a
			d2 = ( runoffij - tau ) * ( 1 - a )
		end	
		d1_directionij = d1_direction[i,j]
		d2_directionij = d2_direction[i,j]
		if d1_directionij == 1 # N
		  	runoff[ i - 2, j ] += d1 
		elseif d1_directionij == 2 #NE
			runoff[ i - 1, j + 1 ] += d1
		elseif d1_directionij == 3 #SE
			runoff[ i + 1, j + 1 ] += d1
		elseif d1_directionij == 4 #S
			runoff[ i + 2, j ] += d1 
		elseif d1_directionij == 5 #SW
			runoff[ i + 1, j - 1 ] += d1
		elseif d1_directionij == 6 #NW
			runoff[ i - 1, j - 1 ] += d1
		end
		if d2_directionij == 1 # N
		  	runoff[ i - 2, j ] += d2 
		elseif d2_directionij == 2 #NE
			runoff[ i - 1, j + 1 ] += d2
		elseif d2_directionij == 3 #SE
			runoff[ i + 1, j + 1 ] += d2
		elseif d2_directionij == 4 #S
			runoff[ i + 2, j ] += d2
		elseif d2_directionij == 5 #SW
			runoff[ i + 1, j - 1 ] += d2
		elseif d2_directionij == 6 #NW
			runoff[ i - 1, j - 1 ] += d2
		end
	end	
return runoff
end



function bucketlocation(bucket, incidence)
	temp = copy(bucket) 
	marker_bad = Array(Int64, 0)
	marker_good = Array(Int64, 0) #  rows (y) 
	for i = 1:34
		if incidence[round(Int64,floor(temp[i,1])), round(Int64,floor(temp[i,2]))] == 0 # check if in the checker
			push!(marker_bad, i) # wierd bucket ID recorded
		else
			push!(marker_good, i)
			temp[i,1] = map(Int64, floor(temp[i,1]))
			temp[i,2] = map(Int64, floor(temp[i,2]))
		end
	end
	good = temp[collect(marker_good),:]
neib_loc = zeros(length(marker_bad), 8)
bad = zeros(length(marker_bad), 2)
marked_bucket = 0
	for j in marker_bad
		marked_bucket += 1
		m = floor(temp[j,1]) 
		n = floor(temp[j,2])
		neib_loc[marked_bucket, 1:2] = [ m - 1, n ] #top
		neib_loc[marked_bucket, 3:4] = [ m, n + 1 ]	# right
		neib_loc[marked_bucket, 5:6] = [ m + 1, n ] # bottom
		neib_loc[marked_bucket, 7:8] = [ m, n - 1 ] #left
		bad[marked_bucket,:] = copy(temp[j,:])
	end
return Array{Int64}[good, bad, neib_loc, marker_good, marker_bad]
end


function bucket_water(runoff_sim_output, bucket_loc, method = "WAVG")
good_buckets = bucket_loc[1]
marker_good = bucket_loc[4]
marker_bad = bucket_loc[5]
water_amount = zeros(34)
iteration = 0
	for i in marker_good
		iteration += 1
		water_amount[i] = runoff_sim_output[good_buckets[iteration,1], good_buckets[iteration,2]]
	end
bad_buckets = bucket_loc[2]
bad_neib_location = bucket_loc[3]
iteration = 0
	if method == "WAVG"
		for j in marker_bad
			iteration += 1
			weight = zeros(4) #[top, bottom, left, right]
			weight[1] = bad_buckets[iteration,1] % 1 # top
			weight[2] = 1 - weight[1] # bot
			weight[3] = bad_buckets[iteration,2] % 1 # left
			weight[4] = 1 - weight[3]  # right
			weight = weight/sum(weight)
			water_amount[j] = runoff_sim_output[bad_neib_location[iteration,1], bad_neib_location[iteration,2]]*weight[1] + #top
						  	runoff_sim_output[bad_neib_location[iteration,3], bad_neib_location[iteration,4]]*weight[4] + #right
						  	runoff_sim_output[bad_neib_location[iteration,5], bad_neib_location[iteration,6]]*weight[2] + # bottom
						  	runoff_sim_output[bad_neib_location[iteration,7], bad_neib_location[iteration,8]]*weight[3] #left
		end
	elseif method == "NNB" 
		for j in marker_bad
			iteration += 1
			i_dec = bad_buckets[iteration, 1] % 1
			j_dec = bad_buckets[iteration, 2] % 1
			if i_dec >= 0.5 && j_dec >= 0.5 # BR
				if i_dec > j_dec # go down
					i_coord = ceil(bad_buckets[iteration, 1])
					j_coord = floor(bad_buckets[iteration, 2])  
				else # go right
					i_coord = floor(bad_buckets[iteration, 1])
					j_coord = ceil(bad_buckets[iteration, 2])
				end
			elseif i_dec < 0.5 && j_dec >= 0.5 # TR
				if (0.5 - i_dec ) > j_dec # go up
					i_coord = floor(bad_buckets[iteration, 1]) - 1
					j_coord = floor(bad_buckets[iteration, 2])
				else # go right
					i_coord = floor(bad_buckets[iteration, 1]) 
					j_coord = ceil(bad_buckets[iteration, 2])
				end
			elseif i_dec < 0.5 && j_dec < 0.5 # TL
				if (0.5 - i_dec) > (0.5 - j_dec) # go up
					i_coord = floor(bad_buckets[iteration, 1]) - 1
					j_coord = floor(bad_buckets[iteration, 2])
				else # go left
					i_coord = floor(bad_buckets[iteration, 1]) 
					j_coord = floor(bad_buckets[iteration, 2]) - 1
				end
			elseif i_dec > 0.5 && j_dec < 0.5 # BR
				if (1 - i_dec) > j_dec # go left
					i_coord = floor(bad_buckets[iteration, 1]) 
					j_coord = floor(bad_buckets[iteration, 2]) - 1
				else # go down
					i_coord = ceil(bad_buckets[iteration, 1])
					j_coord = floor(bad_buckets[iteration, 2])
				end
			end
			water_amount[j] = runoff_sim_output[i_coord, j_coord]	
		end
	end
	return water_amount
end



function lookup_table_sim(tau, a, method, incidence, runoff_directions, runoff_ranking, infiltration, rainfall_range)
m = size(bucket, 1)
n = length(rainfall_range)
LUT = zeros(m,n)#lookup table
col = 0
for i in rainfall_range
	col += 1
	runoff_sim_output = runoff_sim2(runoff_ranking, i, infiltration, incidence, runoff_directions, tau, a)
	tip = bucket_water(runoff_sim_output, bucket_loc, method)
	LUT[:, col] = tip
end
return LUT
end


function missingdata(final_fantasy)
final_fantasy_II = copy(final_fantasy)
final_fantasy_II[25,1:5] = 0
final_fantasy_II[13,19] = 0
final_fantasy_II[31,19] = 0
return final_fantasy_II
end


function driver2(theta, method, incidence, train, rainfall_table, rainfall_range, heightdata)
final_fantasy_IV = Vector(0)
@time runoff_directions = direction_simulation2(heightdata, incidence, theta[1]);
@time infiltration = infiltration_gen(heightdata, theta[3])
@time lookup_table = lookup_table_sim(theta[2], theta[4], method, incidence, runoff_directions, ranked ,infiltration, rainfall_range)
@time final_fantasy = lookup_table * rainfall_table * 1e-6 ######## m3 to L
@time final_fantasy_II = missingdata(final_fantasy)
for period = 1:size(train,2)
	push!(final_fantasy_IV, ecdf(final_fantasy_II[:,period]))
end
return final_fantasy_IV
end 

function observed_ecdf()
final_fantasy_V = Vector(size(train,2))
for i = 1:size(train,2)
	final_fantasy_V[i] = ecdf(train[:,i])
end
return final_fantasy_V
end 
