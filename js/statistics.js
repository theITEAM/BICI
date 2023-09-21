// Functions qhich calculate statistical quantities

"use strict";

/// Given a vector of values this returns statistical information
function get_statistic(vect)
{
	let vec = copy(vect);
	
	let stat = {};

	let n = vec.length;
	if(n == 0){
		stat.mean = "NA";
		stat.CImin = "NA";
		stat.CImax = "NA";
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

	return stat;
}


/// Gets the minimum and maximum of a list
function get_range(vec)
{
	let min = LARGE, max = -LARGE;
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
	let av = 0.0, av2 = 0.0;
	let N = vec.length;
	for(let s = 0; s < N; s++){
		let val = vec[s]; av += val; av2 += val*val;
	}
	
	let num = av2/N - (av/N)*(av/N);
	
	if(num < TINY*TINY || N < 2) return "-";
	
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
	

/// Returns the gelman rubin statistic
function get_Gelman_Rubin_statistic(cha)
{
	let C = cha.length;
	
	let N = LARGE;
	for(let ch = 0; ch < C; ch++){
		if(cha[ch].length < N) N = cha[ch].length;
	}		
		
	if(N == 0){ error =("no data"); return;}
		
	let mu=[], vari=[];
				
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
	return Math.sqrt(((1-1.0/N)*W + B/N)/W);
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
