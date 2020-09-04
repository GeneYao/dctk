#ifndef _RECEIVER_TABLE_HPP
#define _RECEIVER_TABLE_HPP

#include "dctk.hpp"
#include "CellLib.hpp"
#include "Circuit.hpp"
#include <parser-spef.hpp>
#include "si2dr_liberty.h"


class ReceiverTable {
    const liberty_value_data* cap_table;
	int t_size;
	int c_size;
	
	public:
    ReceiverTable(){
        cap_table = nullptr;
        t_size = 0;
        c_size = 0;
    }
	ReceiverTable(const liberty_value_data* receiver_table){
        load_data(receiver_table);
	}
    
    /// not acctually loading data
    void load_data(const liberty_value_data* receiver_table){
		cap_table = receiver_table;
		t_size = receiver_table->dim_sizes[0];
		c_size = receiver_table->dim_sizes[1];
        
    } 
    
	float get_cap(int idx_t, int idx_c){
		if(idx_t<0 || idx_t>t_size || idx_c<0 || idx_c>c_size) return -1;
		return cap_table->values[idx_t*c_size + idx_c];
	}

	float tr_diff(float tr, int idx){
		return tr - cap_table->index_info[0][2*idx-1];
	}

	std::size_t get_tsize(){
		return t_size;
	}
	std::size_t get_csize(){
		return c_size;
	}

	float get_c(float tran, float ceff);
	float get_c1(float ceff);
	float get_c2(float ceff);
	float get_c3(float ceff);
};
#endif