"use strict";
// Functions which calculate statistical quantities

/// Given a vector of values this returns statistical information
function get_statistic(vect)
{
	let vec = copy(vect);

	let n = vec.length;
	
	let stat = {n:n};
	if(n == 0){
		stat.mean = "NA";
		stat.CImin = "NA";
		stat.CImax = "NA";
	} 
	else{
		let min = LARGE, max = -LARGE;
		for(let i = 0; i < vec.length; i++) {
			let val = vec[i];
			if(val < min) min = val;
			if(val > max) max = val;
		}
		
		if(min == max){
			stat.mean = min;
			stat.CImin = min;
			stat.CImax = min;
		}
		else{		
			let sum = 0, sum2 = 0;
			for(let i = 0; i < vec.length; i++) {
				sum += vec[i];
				sum2 += vec[i]*vec[i];
			}
			sum /= n;
			sum2 /= n;

			stat.mean = sum;

			vec.sort(function(a, b){return a - b});

			if(n >= 2) {
				let i = Math.floor((n-1)*0.025);
				let f = (n-1)*0.025 - i;
				stat.CImin = vec[i]*(1-f) + vec[i+1]*f;

				i = Math.floor((n-1)*0.975);
				f = (n-1)*0.975 - i;
				stat.CImax = vec[i]*(1-f) + vec[i+1]*f;
			} 
			else{
				stat.CImin = vec[0];
				stat.CImax = vec[0];
			}
		}
	}

	return stat;
}


/// Gets the minimum and maximum of a list
function get_range(vec)
{
	let min = VLARGE, max = -VLARGE;
	for(let i = 0; i < vec.length; i++){
		let val = vec[i];
		if(val < min) min = val;
		if(val > max) max = val;
	}
	
	return { min:min, max:max};
}


/// Calculate the effective sample size for a list of numbers
function get_effective_sample_size(vec)
{
	let N = vec.length;
	
	let min = LARGE, max = -LARGE;
	for(let i = 0; i < N; i++) {
		let val = vec[i];
		if(val < min) min = val;
		if(val > max) max = val;
	}
	if(min == max || N < 2) return "-";
		
	let av = 0.0, av2 = 0.0;
	for(let s = 0; s < N; s++){
		let val = vec[s]; av += val; av2 += val*val;
	}
	let num = av2/N - (av/N)*(av/N);
	let mean = av/N, sd = Math.sqrt(num);
	
	for(let s = 0; s < N; s++) vec[s] = (vec[s]-mean)/sd;
	
	let sum = 1.0;
	for(let d = 1; d < N/2; d++){
		let a = 0.0; for(let s = 0; s < N-d; s++) a += vec[s]*vec[s+d]; 
		let cor = a/(N-d); if(cor < 0) break;
		sum += 0.5*cor;			
	}
	
	return Math.floor(N/sum);
}	
	

/// Returns the Gelman-Rubin statistic
function get_Gelman_Rubin_statistic(cha)
{
	let C = cha.length;
	
	let N = LARGE;
	for(let ch = 0; ch < C; ch++){
		if(cha[ch].length < N) N = cha[ch].length;
	}		
		
	if(N == 0){ error =("no data"); return;}
		
	let mu=[], vari=[];
	
	{ // Checks to see if values along chain are all equal
		for(let ch = 0; ch < C; ch++){ 
			let val_basic = cha[ch][0];
			let i = 0; while(i < N && cha[ch][i] == val_basic) i++;
			if(i == N) return "-";
		}
	}
	
	let muav = 0.0;
	for(let ch = 0; ch < C; ch++){ 
		let valav = 0.0; for(let i = 0; i < N; i++) valav += cha[ch][i]/N;
		let varr = 0.0; for(let i = 0; i < N; i++) varr += (cha[ch][i]-valav)*(cha[ch][i]-valav)/(N-1);
		mu[ch] = valav;
		vari[ch] = varr;
		muav += mu[ch]/C;
	}
	
	let W = 0.0; for(let ch = 0; ch < C; ch++) W += vari[ch]/C;
	
	let B = 0.0; for(let ch = 0; ch < C; ch++) B += (mu[ch]-muav)*(mu[ch]-muav)*N/(C-1);
	return Math.sqrt(((1-1.0/N)*W + B/N)/W).toPrecision(pre);
}


/// Works out the mean and credible interval for a series of lines
function get_line_stats(line)
{
	if(line.length == 0) error("Problem with line stat");
	
	let T = line[0].length;

	let line_stats = [];
	for(let t = 0; t < T; t++){
		let vec = [];
		for(let i = 0; i < line.length; i++) vec.push(line[i][t].y);
		let stat = get_statistic(vec);
		
		line_stats.push({x:line[0][t].x, y:stat.mean, CImin:stat.CImin, CImax:stat.CImax});
	}

	return line_stats;
}			


/// Works out the mean for a series of lines
function get_line_mean(line)
{
	if(line.length == 0) error("Problem with line stat");
	
	let T = line[0].length;

	let line_stats = [];
	
	let av = [];
	for(let t = 0; t < T; t++) av[t] = 0;
	
	for(let i = 0; i < line.length; i++){
		let li = line[i];
		for(let t = 0; t < T; t++) av[t] += li[t].y;
	}
	
	let num = line.length;
	for(let t = 0; t < T; t++){
		//line_stats.push({x:line[0][t].x, y:av[t]/num});
		line_stats.push(av[t]/num);
	}

	return line_stats;
}			


/// Gets the Pearson correlation coefficient between two sets of measurements
function get_correlation(vecA,vecB)
{
	if(vecA.length != vecB.length){ error("Cannot get correlation"); return 0;}
	
	let avA = 0, avA2 = 0, avB = 0, avB2 = 0, avAB = 0;
	let N = vecA.length;
	for(let i = 0; i < N; i++){
		let valA = vecA[i], valB = vecB[i];
		avA += valA; avA2 += valA*valA;
		avB += valB; avB2 += valB*valB;
		avAB += valA*valB;
	}
	
	let varA = avA2/N - (avA/N)*(avA/N);
	let varB = avB2/N - (avB/N)*(avB/N);
	
	if(varA < TINY || varB < TINY) return 0;
	
	return (avAB/N - (avA/N)*(avB/N))/Math.sqrt(varA*varB);
}


/// Calculates the xi sq distribution - used for testing for in trans (p-val)
function calc_xi_sq()
{
	let dx = 0.1;
	let k = H_BIN-1;
	
	let prob=[];
	let imax = 10*k/dx;
	let sum = 0;
	for(let i = 0; i < imax; i++){
		let x = i*dx;
		prob[i] = Math.exp((k/2-1)*Math.log(x)-x/2);	
		sum += prob[i];
	}
	
	let prob_above=[];
	let sum2 = 0;
	for(let i = imax-1; i >= 0; i--){
		prob_above[i] = sum2;
		sum2 += prob[i]/sum;
	}
	
	return { prob_above:prob_above, dx:dx};
}


/// Gets the p-value from xis
function xi_sq_p_value(xis)
{
	let i = Math.floor(xis/xi_sq.dx);
	if(i >= xi_sq.prob_above.length) return TINY;
	return xi_sq.prob_above[i];
}
