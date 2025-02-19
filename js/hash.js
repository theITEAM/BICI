"use strict";
// Functions which create a hash table

class Hash
{
	store=[];
	num=[];
	list=[];
	ran=[];
	ran2=[];
	ran_base;
	
	create_ran()
	{
		for(let i = 0; i <= TI_DIV_MAX; i++){
			this.ran[i] = Math.random();
			this.ran2[i] = Math.random();
		}
		this.ran_base = Math.random();
	}
	
	
	/// Creates a hash number based on two vector
	create_H(v,n)
	{
		let ba = this.ran_base;
		let ran = this.ran;
		let ran2 = this.ran2;
		
		let H = 0;
		for(let k = 0; k < v.length; k++){
			H += ran[k]*(v[k]+ba) + ran2[k]*n[k];
		}
		return String(H);
	}
	
	/// Given a references finds an element
	find(ref)
	{
		return hash_find(this.store,ref);
	}
	
	
	/// Adds a reference to the hash table
	add(ref,num)
	{
		hash_add(this.store,ref,num);
	}
	
	
	/// Removes a reference to th hash table
	remove(ref)
	{
		hash_remove(this.store,ref);
	}
	
	
	/// These are used to generate phylogenetic 
	set_phylo_num(N)
	{
		for(let i = 0; i < N; i++){
			this.num[i] = Math.round(1000000000000*Math.random());
		}
	}
	
	
	/// Gets a hash code from a phylogenetic vector
	get_phylo_code(vec)
	{
		let nu = 0;
		for(let i = 0; i < vec.length; i++){
			if(vec[i] >= this.num.length){
				pr(vec[i]+" "+this.num.length+" pp");
				error("Out of range hash");
			}
			nu = nu ^ this.num[vec[i]]; 
		}
		return nu;
	}
	
	
	/// Adds an object with a counter to the hash table
	add_phylo_ori(ref,t_vec,mnum_vec,arr)
	{
		let code = hash_get_code(ref);
	
		if(this.store[code] == undefined){
			this.list.push(code);
			this.store[code]=[];
		}
		
		for(let k = 0; k < this.store[code].length; k++){
			if(this.store[code][k].ref == ref){
				this.store[code][k].num++;
				for(let ii = 0; ii < t_vec.length; ii++){
					this.store[code][k].t_sum[ii] += t_vec[ii];
					this.store[code][k].mnum_sum[ii] += mnum_vec[ii];
				}
				return;
			}
		}
		
		this.store[code].push({ref:ref, t_sum:t_vec, mnum_sum:mnum_vec, arr:arr, num:1});
	}

	
	/// Returns the most frequent object in the table
	get_phylo_ori()
	{
		let num_tot = 0;
		let num_max = 0;
		let st_max;
		for(let k = 0; k < this.list.length; k++){
			let st = this.store[this.list[k]];
			for(let j = 0; j < st.length; j++){
				let stt = st[j];	
				if(stt.num > num_max){ num_max = stt.num; st_max = stt;}
				num_tot += stt.num;
			}
		}
		
		if(num_tot == 0) return;
		
		let t_vec = st_max.t_sum;
		let mnum_vec = st_max.mnum_sum;
		for(let i = 0; i < t_vec.length; i++){
			t_vec[i] /= st_max.num;
			mnum_vec[i] /= st_max.num;
		}
		
		return {freq:st_max.num/num_tot, t_vec:t_vec, mnum_vec:mnum_vec, arr:st_max.arr};
	}
	
	
	/// Adds an object with a counter to the hash table
	add_phylo_br(ref,ob)
	{
		let code = hash_get_code(ref);
	
		if(this.store[code] == undefined){
			this.list.push(code);
			this.store[code]=[];
		}
		
		let k;
		for(k = 0; k < this.store[code].length; k++){
			if(this.store[code][k].ref == ref){
				break;
			}
		}
		
		if(k == this.store[code].length){
			this.store[code].push({ref:ref, pos:[]});
		}
		
		let pos = this.store[code][k].pos;
	
		switch(ob.ty){
		case "OBS":
			{
				let j = 0; while(j < pos.length && !(pos[j].ty == "OBS" && pos[j].m == ob.m)) j++;
				
				if(j < pos.length){
					let po = pos[j];
					po.num++;
					po.t_sum += ob.t;
					po.mnum_sum += ob.mnum;
					let ii = 0; 
					while(ii < po.ind_list.length && po.ind_list[ii].ind != ob.ind) ii++;
					if(ii < po.ind_list.length) po.ind_list[ii].num++;
					else po.ind_list.push({num:1,ind:ob.ind});
				}
				else{
					let add = {ty:"OBS", m:ob.m, num:1, t_sum:ob.t, mnum_sum:ob.mnum, ind_list:[]};
					add.ind_list.push({num:1,ind:ob.ind});
					pos.push(add);
				}
			}
			break;
			
		case "SPLIT":
			{
				let j = 0;
				while(j < pos.length && !(pos[j].ty == "SPLIT" && pos[j].code_br == ob.code_br)) j++;
				
				if(j < pos.length){
					let po = pos[j];
					po.num++;
					po.t_sum += ob.t;
					po.mnum_sum += ob.mnum;
					let ii = 0; 
					while(ii < po.ind_list.length && po.ind_list[ii].ind != ob.ind) ii++;
					if(ii < po.ind_list.length) po.ind_list[ii].num++;
					else po.ind_list.push({num:1,ind:ob.ind});
				}
				else{
					let add = {ty:"SPLIT", code_br:ob.code_br, list_br:ob.list_br, list_br2:ob.list_br2, num:1, t_sum:ob.t, mnum_sum:ob.mnum, ind_list:[]};
					add.ind_list.push({num:1,ind:ob.ind});
					pos.push(add);
				}
			}
			break;	
		}
	}
	
	
	/// Gets the next phylogenetic branch
	get_phylo_br(obs_list,t_now)
	{
		let code_tot = this.get_phylo_code(obs_list);
		let ref = String(code_tot);
		
		let code = hash_get_code(ref);
		
		if(this.store[code] == undefined) return;
		
		let k;
		for(k = 0; k < this.store[code].length; k++){
			if(this.store[code][k].ref == ref){
				break;
			}
		}
		
		if(k == this.store[code].length) return; 
		
		let pos = this.store[code][k].pos;
	
		let max = 0;
		let pos_sel;
		let num_tot = 0;
		for(let j = 0; j < pos.length; j++){
			let po = pos[j];
			if(po.t_sum/po.num > t_now){
				if(po.num > max){ max = po.num; pos_sel = po;}
				num_tot += pos[j].num;
			}
		}
		
		if(pos_sel){
			pos_sel.frac = pos_sel.num/num_tot;
		}
		
		return pos_sel;
	}
	
	
	/// Prints the possible values for splitting
	print_pos()
	{
		pr(this.list.length+" len");
		for(let i = 0; i < this.list.length; i++){
			let st = this.store[this.list[i]];
			for(let j = 0; j < st.length; j++){
				let stt = st[j];
				pr("SPLIT");
				pr(stt);
				for(let k = 0; k < stt.pos.length; k++){
					let po = stt.pos[k];
					pr("pos");
					pr(po);
				}
			}
		}
	}
}

// Hash functions used to find compartments and transitions

/// Converts from a reference to a code which is between zero and HASH_MAX
function hash_get_code(ref)
{
	let hash = 0;

	if(ref.length == 0) return hash;

	for(let i = 0; i < ref.length; i++){
		let ch = ref.charCodeAt(i);
		hash = ((hash << 5) - hash) + ch;
		hash = hash & hash;
	}

	if(hash < 0) hash = -hash;
	
	return hash%HASH_MAX;
}
	
	
/// Given a references finds an element
function hash_find(store,ref)
{
	let code = hash_get_code(ref);
	if(store[code] != undefined){
		let sto = store[code];
		for(let i = 0; i < sto.length; i++){
			if(sto[i].ref == ref) return sto[i].num;
		}
	}
}
	
	
/// Adds a reference to the hash table
function hash_add(store,ref,num)
{
	let code = hash_get_code(ref);
	if(store[code] == undefined) store[code]=[];
	
	let vec = store[code];
	for(let k = 0; k < vec.length; k++){
		if(vec[k].ref == ref){
			if(vec[k].num != num) error("hash num is wrong");
			return;
		}
	}
	
	store[code].push({ref:ref,num:num});
}
	
	
/// Removes a reference to the hash table
function hash_remove(store,ref)
{
	let code = hash_get_code(ref);
	if(store[code] == undefined) store[code]=[];
	
	let vec = store[code];
	for(let k = 0; k < vec.length; k++){
		if(vec[k].ref == ref){
			vec.splice(k,1);
			return;
		}
	}
}


/// Redos hash table
function hash_redo(hash,comp)
{
	hash.length = 0;
	for(let c = 0; c < comp.length; c++){
		hash_add(hash,comp[c].name,c);
	}
}
