
#include "dctk.hpp"
#include "CellLib.hpp"
#include "Circuit.hpp"
#include "receiver_table.hpp"
#include <parser-spef.hpp>
#include "si2dr_liberty.h"

float ReceiverTable::get_c(float tran, float ceff){
	int idx_t=-1, idx_c=-1;
	
	for(int i=0; i<c_size-1; ++i){
		if( (tran >= cap_table->index_info[0][i]) && (tran <= cap_table->index_info[0][i+1]) ) idx_t=i;
	}
	
	if(idx_t==-1 && tran < cap_table->index_info[0][0]) idx_t=0;
	if(idx_t==-1 && tran > cap_table->index_info[0][c_size - 2]) idx_t=c_size-2;
	
	for(int i=0; i<c_size-1; ++i){
		if( (ceff >= cap_table->index_info[1][i]) && (ceff <= cap_table->index_info[1][i+1]) ) idx_c=i;
	}
	if(idx_c==-1 && ceff < cap_table->index_info[1][0]) idx_c=0;
	if(idx_c==-1 && ceff > cap_table->index_info[1][c_size-2]) idx_c=c_size-2;
	
		std::cout << "receiver " << idx_t << " " << idx_c << std::endl;
		
	float c1 = get_cap(idx_t,idx_c) + ( get_cap(idx_t+1,idx_c) - get_cap(idx_t,idx_c) ) * ( tran -cap_table->index_info[0][idx_t] ) / ( cap_table->index_info[0][idx_t+1] - cap_table->index_info[0][idx_t] );
	float c2 = get_cap(idx_t,idx_c+1) + ( get_cap(idx_t+1,idx_c+1) - get_cap(idx_t,idx_c+1) ) * ( tran -cap_table->index_info[0][idx_t] ) / ( cap_table->index_info[0][idx_t+1] - cap_table->index_info[0][idx_t] );
	float cc = c1 + (c2-c1) * ( ceff-cap_table->index_info[1][idx_c] ) / ( cap_table->index_info[1][idx_c+1] - cap_table->index_info[1][idx_c] );
		std::cout << "c1 " << c1 << " c2 " << c2 << " cc " << cc << std::endl;
	return cc;
}
float ReceiverTable::get_c1(float ceff){
	
	int idx_t=1, idx_c=-1;
	
	for(int i=0; i<6; ++i){
		if( (ceff >= cap_table->index_info[1][i]) && (ceff <= cap_table->index_info[1][i+1]) ) idx_c=i;
	}
	if(idx_c==-1 && ceff < cap_table->index_info[1][0]) idx_c=0;
	if(idx_c==-1 && ceff > cap_table->index_info[1][c_size-2]) idx_c=c_size-2;
		
	return get_cap(idx_t,idx_c) + ( get_cap(idx_t,idx_c+1) - get_cap(idx_t,idx_c) ) * ( ceff -cap_table->index_info[1][idx_c] ) / ( cap_table->index_info[1][idx_c+1] - cap_table->index_info[1][idx_c] );
		
}

float ReceiverTable::get_c2(float ceff){
	
	int idx_t=3, idx_c=-1;
	
	for(int i=0; i<c_size-1; ++i){
		if( (ceff >= cap_table->index_info[1][i]) && (ceff <= cap_table->index_info[1][i+1]) ) idx_c=i;
	}
	if(idx_c==-1 && ceff < cap_table->index_info[1][0]) idx_c=0;
	if(idx_c==-1 && ceff > cap_table->index_info[1][c_size-2]) idx_c=c_size-2;
		
	return get_cap(idx_t,idx_c) + ( get_cap(idx_t,idx_c+1) - get_cap(idx_t,idx_c) ) * ( ceff -cap_table->index_info[1][idx_c] ) / ( cap_table->index_info[1][idx_c+1] - cap_table->index_info[1][idx_c] );
		
	
}

float ReceiverTable::get_c3(float ceff){
	
	int idx_t=5, idx_c=-1;
	
	for(int i=0; i<c_size-1; ++i){
		if( (ceff >= cap_table->index_info[1][i]) && (ceff <= cap_table->index_info[1][i+1]) ) idx_c=i;
	}
	if(idx_c==-1 && ceff < cap_table->index_info[1][0]) idx_c=0;
	if(idx_c==-1 && ceff > cap_table->index_info[1][c_size-2]) idx_c=c_size-2;
		
	return get_cap(idx_t,idx_c) + ( get_cap(idx_t,idx_c+1) - get_cap(idx_t,idx_c) ) * ( ceff -cap_table->index_info[1][idx_c] ) / ( cap_table->index_info[1][idx_c+1] - cap_table->index_info[1][idx_c] );
		
	
}
