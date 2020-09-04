#include <fstream>
#include <algorithm>
#include <iostream>
#include "current_table.hpp"

#include "delay_calculator.hpp"
#include "dctk.hpp"
#include "CellLib.hpp"
#include "Circuit.hpp"
#include <parser-spef.hpp>
#include "si2dr_liberty.h"


void CurrentTableData::read_table(const std::vector<liberty_value_data*>& table_in)
{
    
    int len = table_in.size();
    for(int i = 0; i < len; i++) {            
        liberty_value_data* table;
        table = table_in[i];
		
        int num_dimensions = table -> dimensions;
        int total_element_count = 1;
        int ceff_index = i % ceff_index_size;
        int tran_index = i / tran_index_size;
        for (int d = 0; d < num_dimensions; d++) {
            // index values 
            int dim_size = table->dim_sizes[d];
            LONG_DOUBLE *index_values = table->index_info[d];
            if (d == 0) {
                tran[tran_index] = index_values[0];
            }
            else if (d == 1) {
                ceff[ceff_index] = index_values[0];
            }
            else if (d == 2) {
                for (int i = 0; i<dim_size; i++) {
                    time[tran_index][ceff_index].push_back(index_values[i]);
				}	
            }
            total_element_count = total_element_count * dim_size;
        }
		
        // array of values 
        for (int i = 0; i < total_element_count; i++) {
            iout[tran_index][ceff_index].push_back(table->values[i]);
        }
		
    }
}


void CurrentTable::Search(float input_t, float input_c) {
    _input_t = input_t;
    _input_c = input_c;
    // --- search specific tables ---
    int t = -1, c = -1;
	

    for (int i = 0; i <= tran_index_size - 2; i++) {
		
        if (input_t >= tran[i] && input_t <= tran[i+1]) {
            t = i;
            break;
        }       
		
    }
    for (int i = 0; i <= ceff_index_size - 2; i++) {
		
        if (input_c >= ceff[i] && input_c <= ceff[i+1]) {
            c = i;
            break;
        }        
		
    }
	
    if (c == -1) {
		
		if( input_c > ceff[ceff_index_size - 1] ){
			c = ceff_index_size - 2;
		}else{
			c = 0;
		}
    }
    target_table_index = std::make_pair(t, c);
	
}

void CurrentTable::interpolate_1D(const std::vector<float>& y,const std::vector<float>& t
				   ,std::vector<float>& qy, std::vector<float>& qt
				   ,float dt)
{
	/*
		linear interpolate to points (t(i), y(i))
		assume t is sorted in increasing order...
		y: input y(i)
		t: input t(i)
		qy: interpolated y(j)
		qt: interpolated t(j) from 0 to max(t), step = dt
		dt: step size of interpolated t(j)
	*/
	
	
	int tsize = t.size();
    
	int qsize = ceil( t.back() /dt );
    
	qy.resize(qsize);
	qt.resize(qsize);
	int i=0;
	for(int j=0; j<qsize; ++j){
		if(j>0) qt[j] = qt[j-1]+dt;
		else    qt[j] = 0;
		
		if(qt[j]<t[0]){
			// before t[0]
			qy[j] = 0;
		}else if(i>=tsize-1){
			// after t[back]
			qy[j]=0;
		}else{
			while( i<tsize-1 && qt[j]>=t[i+1] ){
				i=i+1; // move to next interval
			}
			qy[j] = y[i] + (qt[j]-t[i]) / (t[i+1]-t[i]) * (y[i+1]-y[i]);
		}
		
		
	}
}
 
void CurrentTable::waveform_projection(Waveform* i_composite){
	/// find composite waveform using projection matrices
	
	/// combine 2 waveforms with the same Ceff but different Tr
	
	int t1_index = target_table_index.first, c1_index = target_table_index.second;
	get_waveform(i_t1_c1,t1_index,c1_index);
	get_waveform(i_t2_c1,t1_index+1,c1_index);
	get_waveform(i_t1_c2,t1_index,c1_index+1);
	get_waveform(i_t2_c2,t1_index+1,c1_index+1);
	i_t1_c1.trace_current(project_tail_point);
	i_t2_c1.trace_current(project_tail_point);
	i_t1_c2.trace_current(project_tail_point);
	i_t2_c2.trace_current(project_tail_point);
	
	float c1 = ceff[c1_index], c2 = ceff[c1_index+1];
	i_tr_c1.set_ceff(c1);
	i_tr_c2.set_ceff(c2);
	i_tr_c1.set_polarity( this->get_polarity() );
	i_tr_c2.set_polarity( this->get_polarity() );
	
	
	/// Binary search for Tran by gate delay
	
	set_tr(0.5);
	set_tr_max(1.0);
	set_tr_min(0.0);
	float td1 = i_t1_c1.get_tdelay();
	float td2 = i_t2_c1.get_tdelay();
	float tran1 = tran[t1_index];
	float tran2 = tran[t1_index+1];
	float tran_ratio = (_input_t-tran1)/(tran2-tran1);
	float td_target = td1 + (td2-td1)*tran_ratio;
	float td;
	
#ifdef DELAY_CAL_DEBUGGING_MODE
	if(tran_ratio>1.0){
		std::cout << "Warning: <WaveProject> Td is large than Lib limit, using extrapolation" << std::endl;
	}
#endif
	
	/// c1
	for(int i=0; i<max_ite_t1; ++i){
		
		i_tr_c1.get_projection(&i_t1_c1, &i_t2_c1, tr);
		td = i_tr_c1.get_tdelay();
				
		if(td_target==td){
			break;
		}else if(td_target<td){
			tr_max = tr;
		}else{
			tr_min = tr;
		}
		
		tr = (tr_max+tr_min)/2;
	}
	
#ifdef DELAY_CAL_DEBUGGING_MODE
	if(std::abs(td-td_target)>0.01){
		std::cout << "@@ Warning: <WaveProject> searching for td1 does not converge" << std::endl;
		std::cout << "@@ target is " << td_target << ", found " << td << std::endl;
	}
#endif
	
	set_tr(0.5);
	set_tr_max(1.0);
	set_tr_min(0.0);
	td1 = i_t1_c2.get_tdelay();
	td2 = i_t2_c2.get_tdelay();
	td_target = td1 + (td2-td1)*tran_ratio;
	
	
	/// c2
	for(int i=0; i<max_ite_t2; ++i){
		
		i_tr_c2.get_projection(&i_t1_c2, &i_t2_c2, tr);	
		td = i_tr_c2.get_tdelay();
		
		if(td_target==td){
			break;
		}else if(td_target<td){
			tr_max = tr;
		}else{
			tr_min = tr;
		}
		
		tr = (tr_max+tr_min)/2;
	}
	
#ifdef DELAY_CAL_DEBUGGING_MODE
	if(std::abs(td-td_target)>0.01){
		std::cout << "@@ Warning: <WaveProject> searching for td2 does not converge" << std::endl;
		std::cout << "@@ target is " << td_target << ", found " << td << std::endl;
	}
#endif
	
	
	/// Binary search for Ceff by total charge
	set_cr(0.5);
	set_cr_max(1.0);
	set_cr_min(0.0);
	float q_target = _input_c * vdd, q;
	
	if(_input_c > get_ceff_max()){
#ifdef DELAY_CAL_DEBUGGING_MODE
		std::cout << "Warning: <WaveProject> Ceff is larger than Lib limit, using extrapolation" << std::endl;
#endif

		set_cr(2.0);
		set_cr_max(3.0);
		set_cr_min(1.0);
	}
	
	if(_input_c < get_ceff_min()){
#ifdef DELAY_CAL_DEBUGGING_MODE
		std::cout << "Warning: <WaveProject> Ceff is smaller than Lib limit, using extrapolation" << std::endl;
#endif

		set_cr(-1.0);
		set_cr_max(0.0);
		set_cr_min(-2.0);
	}


	for(int i=0; i<max_ite_cap; ++i){
		
		i_composite->get_projection(&i_tr_c1,&i_tr_c2,cr);
		q = std::abs( i_composite->total_charge() );
				
		if(q_target==q){
			break;
		}else if(q_target<q){
			cr_max = cr;
		}else{
			cr_min = cr;
		}
		
		cr = (cr_max+cr_min)/2;
	}
	
#ifdef DELAY_CAL_DEBUGGING_MODE
	if(std::abs(q-q_target)>0.01){
		std::cout << "@@ Warning: <WaveProject> searching for charge does not converge" << std::endl;
		std::cout << "@@ target is " << q_target << ", found " << q << std::endl;
	}

		/// dump waveforms
		i_t1_c1.dump_gnu("i_t1_c1.it");
		i_t2_c1.dump_gnu("i_t2_c1.it");
		i_tr_c1.dump_gnu("i_tr_c1.it");
		i_t1_c2.dump_gnu("i_t1_c2.it");
		i_t2_c2.dump_gnu("i_t2_c2.it");
		i_tr_c2.dump_gnu("i_tr_c2.it");
		i_composite->dump_gnu("i_tr_cc.it");

#endif
	
}

float CurrentTable::get_waveform(Waveform& waveform
				 ,float dt, float input_t, float input_c)
{
	this->Search(input_t, input_c);
	
	if(dt==0){
		dt = ( time[target_table_index.first+1][target_table_index.second + 1].back() )/ this->num_step; 
	}
	/// time-step based waveform		 
	/// this->table_interpolation(waveform.i_data(), waveform.t_data(), dt);	
	
	/// projection based waveform
	this->waveform_projection(&i_composite);
	
	interpolate_1D(i_composite.i(),i_composite.t(),waveform.i_data(),waveform.t_data(),dt);
	
	/// set waveform properties
	waveform.set_ceff(input_c);
	waveform.set_dt(dt);
	waveform.set_polarity( this->get_polarity() );
	
	return dt;
}

void CurrentTable::get_waveform(Waveform& waveform , int i, int j)
{
	/// i,j are Tr/Ceff indexes, return a specific wavefrom in libaray that is not interpolated.
	if(i < 0 || i > ceff_index_size - 1 || j < 0 || j > tran_index_size - 1) return;
    
	waveform.i_data() = iout[i][j];
	waveform.t_data() = time[i][j];
	
	waveform.set_ceff(ceff[j]);
	waveform.set_polarity( this->get_polarity() );
}

float CurrentTable::lookup_td(float td, float tr, float dt, float ratio){
		
	int ite=0;
	
	// use pointer to avoid copy/swap of waveforms
	Waveform* upper_ptr = &i_upper;
	Waveform* lower_ptr = &i_lower;
	Waveform*  find_ptr = &i_find;
	float c_upper = get_ceff_max();
	float c_lower = get_ceff_min();
	
	get_waveform( *upper_ptr, dt, tr, c_upper);
	get_waveform( *lower_ptr, dt, tr, c_lower);
	
	float c_find = std::sqrt(c_upper*c_lower);
	float t_find;
	while(true){
		get_waveform(*find_ptr, dt, tr, c_find);
		t_find = find_ptr->get_tdelay();
		
		if( td == t_find ){
			break;			
		}else if( td < t_find ){
			std::swap(upper_ptr,find_ptr);
			c_upper = c_find;
		}else{
			std::swap(lower_ptr,find_ptr);
			c_lower = c_find;
		}
		
		c_find = std::sqrt(c_upper*c_lower);
		if(ite++ == max_ite_lookup) break;
	}
		
	return c_find;
}
