#include <iostream>
#include <thread>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <fstream>
#include <algorithm>
#include "current_table.hpp"
#include <stdio.h>
#include <string.h>

#include "si2dr_liberty.h"
#include "dctk.hpp"
#include "CellLib.hpp"
#include "Circuit.hpp"
#include "PiModel.hpp"
#include <parser-spef.hpp>
#include <vector>

#include "current_table.hpp"
#include "receiver_table.hpp"
#include "waveform.hpp"
#include "circuit_sim.hpp"
#include "os_service.h"
#include <fstream>
#include <cassert>
#include <cstring>

#include <chrono>
#include "delay_calculator.hpp"
#include <omp.h>


void DelayCalculator::read_config(char* filename){
	
	if( filename==nullptr ){ 
        std::cout << "No config file found, using default configuration." << std::endl;
		return;
	}
    
	std::ifstream ifs(filename);
	if( !ifs ){ 
        std::cout << "No config file found, using default configuration." << std::endl;
		return;
	}
    
    std::cout << "load from configuration file: " << filename << std::endl;
	while( !ifs.eof() ){
		std::string option;
		float value;
		ifs >> option >> value;
		
		if(option=="search_range_ceff1_min"){
			set_cmin_range(value);
			std::cout << "<" << filename << "> search_range_ceff1_min = " << value << std::endl;
		}else if(option=="search_range_ceff1_max"){
			set_cmax_range(value);
			std::cout << "<" << filename << "> search_range_ceff1_max = " << value << std::endl;
		}else if(option=="search_max_iteration_s1"){
			set_max_ite_s1(value);
			std::cout << "<" << filename << "> search_max_iteration_s1 = " << value << std::endl;
		}else if(option=="project_number_step"){
			set_project_number_step(value);
			std::cout << "<" << filename << "> project_number_step = " << value << std::endl;
		}else if(option=="project_max_iteration_t1"){
			set_max_ite_t1(value);
			std::cout << "<" << filename << "> project_max_iteration_t1 = " << value << std::endl;
		}else if(option=="project_max_iteration_t2"){
			set_max_ite_t2(value);
			std::cout << "<" << filename << "> project_max_iteration_t2 = " << value << std::endl;
		}else if(option=="project_max_iteration_cap"){
			set_max_ite_cap(value);
			std::cout << "<" << filename << "> project_max_iteration_cap = " << value << std::endl;
		}else if(option=="project_tail_point"){
			set_project_tail_point(value);
			std::cout << "<" << filename << "> project_tail_point = " << value << std::endl;
		}else if(option=="lookup_td_max_iteration"){
			set_max_ite_lookup(value);
			std::cout << "<" << filename << "> lookup_td_max_iteration = " << value << std::endl;
		}else if(option=="circuit_sim_tolerance"){
			set_sim_tolerance(value);
			std::cout << "<" << filename << "> circuit_sim_tolerance = " << value << std::endl;
		}else if(option=="circuit_sim_step_limit"){
			set_sim_step_limit(value);
			std::cout << "<" << filename << "> circuit_sim_step_limit = " << value << std::endl;
		}else if(option=="cat_ratio_rise"){
			set_cat_ratio_rise(value);
			std::cout << "<" << filename << "> cat_ratio_rise = " << value << std::endl;
		}else if(option=="cat_ratio_fall"){
			set_cat_ratio_fall(value);
			std::cout << "<" << filename << "> cat_ratio_fall = " << value << std::endl;
		}else if(option=="boost_mode"){
			set_boost_mode(value);
			std::cout << "<" << filename << "> boost_mode = " << value << std::endl;
		}else if(option=="num_feature"){
			set_num_feature(value);
			std::cout << "<" << filename << "> num_feature = " << value << std::endl;
		}else if(option=="num_thread"){
			set_num_thread(value);
			std::cout << "<" << filename << "> num_thread = " << value << std::endl;
		}else if(option=="output_feature"){
			set_output_feature(value);
			std::cout << "<" << filename << "> output_feature = " << value << std::endl;
		}
	}
	
	ifs.close();
}


// From https://helloacm.com/cc-function-to-compute-the-bilinear-interpolation/
inline float 
BilinearInterpolation(float q11, float q12, float q21, float q22, float x1, float x2, float y1, float y2, float x, float y) 
{
    const float x2x1 = x2 - x1;
    const float y2y1 = y2 - y1;
    const float x2x = x2 - x;
    const float y2y = y2 - y;
    const float yy1 = y - y1;
    const float xx1 = x - x1;
    return 1.0 / (x2x1 * y2y1) * (
        q11 * x2x * y2y +
        q21 * xx1 * y2y +
        q12 * x2x * yy1 +
        q22 * xx1 * yy1
    );
}

float lookup_table_value(const liberty_value_data* table, float slew, float load) {

    // first check that slew and load are within the table ranges
    unsigned int num_slews = table->dim_sizes[0];
    float min_slew = table->index_info[0][0];
    float max_slew = table->index_info[0][num_slews-1];   

    if ((slew < min_slew) || (slew > max_slew)) {
        std::cout << "Error: slew outside table ranges:  " << std::endl;
        std::cout << "slew = " << slew << "; min_slew = " << min_slew << "; max_slew = " << max_slew << std::endl;
        return 0.0;
    }
    
    unsigned int num_loads = table->dim_sizes[1];
    float min_load = table->index_info[1][0];
    float max_load = table->index_info[1][num_loads-1];   

    if ((load < min_load) || (load > max_load)) {
        std::cout << "Error: load outside table ranges:  " ;
        std::cout << "load = " << load << "; min_load = " << min_load << "; max_load = " << max_load << std::endl;
        return 0.0;
    }


    // x axis == slew
    // y axis == load

    float x2 = 0.0;
    
    unsigned int slew_index = 0;
    for ( ; slew_index < num_slews; slew_index++ ) {
        x2 = table->index_info[0][slew_index];
        if (x2 > slew) {
            break;
        }
    }
    slew_index--;
    float x1 = table->index_info[0][slew_index];
    
    float y2 = 0.0;
    
    unsigned int load_index = 0;
    for ( ; load_index < num_loads; load_index++ ) {
        y2 = table->index_info[1][load_index];
        if (y2 > load) {
            break;
        }
    }
    
    load_index--;
    float y1 = table->index_info[1][load_index];
        
    // now do the table lookup
    unsigned int index = num_slews * slew_index + load_index;
    const float q11 = table->values[index];

    index = num_slews * (slew_index + 1) + load_index;
    const float q21 = table->values[index];

    index = num_slews * slew_index + (load_index + 1);
    const float q12 = table->values[index];

    index = num_slews * (slew_index + 1) + (load_index + 1);
    const float q22 = table->values[index];

    float delay = BilinearInterpolation(q11, q12, q21, q22, x1, x2, y1, y2, slew, load);

    return delay;
}

// this is sort of slow ... how can we make it faster?
spef::Net* find_net(spef::Spef *spef, const std::string& name) {

    size_t numNets = spef->nets.size();
    for (size_t i = 0; i < numNets; i++) {
        if (spef->nets[i].name == name) {
            return &(spef->nets[i]);
        }

    }
    return nullptr;
}


void compute_delays_bilinear(dctk::CellLib* cell_lib, dctk::Circuit* circuit, spef::Spef* spef) {

    // get slew in ps
    const std::string& waveform = circuit->get_input_waveform();
    
    std::istringstream iss(waveform);
    std::vector<std::string> results((std::istream_iterator<std::string>(iss)),
                                     std::istream_iterator<std::string>());

    char *end;
    float slew = strtod(results[1].c_str(), &end);

    // get load from SPEF
    const std::string& netname = circuit->get_name();
    spef::Net* net = find_net(spef, netname);
    float load = net->lcap; // in ff
    
    // driver cell
    const std::string& driver_celltype = circuit->get_driver_celltype();
    // get pointer to driver
    dctk::Cell* cell = cell_lib->get_cell(driver_celltype);

    // get output pin
    dctk::CellPin* output_pin = cell->get_output_pin();

    // get input pin name from driver string
    const std::string& driver = circuit->get_driver() ;
    std::vector<std::string> results1 = split(driver, '/');

    // find arc (second token is input pin)
    dctk::CellArc* arc = output_pin->find_arc(results1[1]);

    // convert slew and load to library units
    float scale_to_ps = cell_lib->get_scale_to_ps();
    slew = slew / scale_to_ps;
    load = load / cell_lib->get_scale_to_ff();

    float rise_delay = lookup_table_value(arc->get_cell_rise_table(), slew, load);
    float fall_delay = lookup_table_value(arc->get_cell_fall_table(), slew, load);
    float rise_transition = lookup_table_value(arc->get_rise_transition_table(), slew, load);
    float fall_transition = lookup_table_value(arc->get_fall_transition_table(), slew, load);

    // save results
    circuit->set_ccs_rise_delay(rise_delay * scale_to_ps);
    circuit->set_ccs_fall_delay(fall_delay * scale_to_ps);
    circuit->set_ccs_rise_slew(rise_transition * scale_to_ps);
    circuit->set_ccs_fall_slew(fall_transition * scale_to_ps);

    return ;
}


void DelayCalculator::load_current_data(dctk::CellLib* cell_lib){
    
    const auto& lib_cell_map = cell_lib->get_cellMap();
    int i = 0;
    
    for( auto ite = lib_cell_map.begin(); ite != lib_cell_map.end(); ite++ ){
        
        cell_map[ite->first] = i;
        
        dctk::Cell* cell = ite->second;
        std::string& cell_name = cell->get_name();
        std::string input_pin_name = "A";      /// kind of not knowing how to get pin name
        
        CurrentTableData* current_rise_ptr = new CurrentTableData( find_current_rise(cell_lib,cell_name,input_pin_name) );
        CurrentTableData* current_fall_ptr = new CurrentTableData( find_current_fall(cell_lib,cell_name,input_pin_name) );
		current_data_rise[i] = current_rise_ptr ;
		current_data_fall[i] = current_fall_ptr ;
        ++i;
    }
}

void DelayCalculator::booster(dctk::CircuitPtrVec* circuitMgr){
    
    if(boost_mode==3){
		
	}else if(boost_mode!=0){
        
        int numCircuits = circuitMgr->size();       
        /// allocate memory for all threads
        float* input_features = new float[num_thread*num_feature];
                    
        #pragma omp parallel 
        {  
        #pragma omp for nowait
        for(int i=0; i<numCircuits; ++i){  
        
            int my_thread = omp_get_thread_num();
            float* input_feature = &( input_features[my_thread*num_feature] );
            std::memset( input_feature, 0.0, num_feature*sizeof(float) );
            
            dctk::Circuit* circuit = (*circuitMgr)[i];
            get_boost_feature(input_feature, circuit);
            
            boost_rise_delay(input_feature, circuit);
            boost_fall_delay(input_feature, circuit);
            boost_rise_slew (input_feature, circuit);
            boost_fall_slew (input_feature, circuit);
        
        } // end of all circuits.
        }
            
        delete [] input_features;
    }
}


void DelayCalculator::get_boost_feature(float* input_feature, dctk::Circuit* circuit ){
    /*
    */
    
    /// get cell type (encoded 0 ~ 24)
    int driver_cell_type = get_encode( circuit->get_driver_celltype() );
    /// get load cell type (encoded 0 ~ 24)
    int load_cell_type = get_encode( circuit->get_load_celltype() );
    /// get RC
    dctk::PiModel interconnect = circuit->get_load_interconnect();
    float cap = interconnect.get_cnear();
    float res= interconnect.get_res();
    /// get input slew
    std::vector<std::string> input_waveform = split( circuit->get_input_waveform(), ' ');
    float tr = std::stof( input_waveform[1] );
    
    input_feature[0] = circuit->get_ccs_rise_delay();
    input_feature[1] = circuit->get_ccs_fall_delay();
    input_feature[2] = circuit->get_ccs_rise_slew();
    input_feature[3] = circuit->get_ccs_fall_slew();
    input_feature[4] = res;
    input_feature[5] = cap;
    input_feature[6] = tr;
    input_feature[7] = res*cap;
    input_feature[8+driver_cell_type] = 1;
    input_feature[33+load_cell_type] = 1;
    
}

void DelayCalculator::boost_rise_delay(float* input_feature, dctk::Circuit* circuit){
                
        float boosted = input_feature[0];
        if(boost_mode==1)
            boosted = score_delay_rise(input_feature);
        else if(boost_mode==2)
            boosted = score_delay_rise_lin(input_feature);
		circuit->set_ccs_rise_delay( boosted );
}

void DelayCalculator::boost_fall_delay(float* input_feature, dctk::Circuit* circuit){
        
        float boosted = input_feature[1];
        if(boost_mode==1)
            boosted = score_delay_fall(input_feature);
        else if(boost_mode==2)
            boosted = score_delay_fall_lin(input_feature);
		circuit->set_ccs_fall_delay( boosted );
}

void DelayCalculator::boost_rise_slew(float* input_feature, dctk::Circuit* circuit){
        
        float boosted = input_feature[2];
        if(boost_mode==1)
            boosted = score_slew_rise(input_feature);
        else if(boost_mode==2)
            boosted = score_slew_rise_lin(input_feature);
		circuit->set_ccs_rise_slew( boosted );
}

void DelayCalculator::boost_fall_slew(float* input_feature, dctk::Circuit* circuit){
        
        float boosted = input_feature[3];
        if(boost_mode==1)
            boosted = score_slew_fall(input_feature);
        else if(boost_mode==2)
            boosted = score_slew_fall_lin(input_feature);
		circuit->set_ccs_fall_slew( boosted );
}

/* Compute delay Stage 1 */
float DelayCalculator::compute_delay_s1(
            CurrentTable& current, ReceiverTable& receiver_1, ReceiverTable& receiver_2, dctk::Circuit* circuit, 
            float& gate_delay, float& intr_delay, float& trans_time, float ratio, int verbose, bool search_delay
){
    /// if search_delay==true
    /// gate_delay is input, golden gate delay
    /// intr_delay is input, original ceff1


#ifdef DELAY_CAL_DEBUGGING_MODE
	if(verbose>=1){
		std::cout << "Stage 1: waveform type = " << ((current.get_polarity()==POS)? "RISE":"FALL")
				  << ", ratio = " << ratio << std::endl;
	}
#endif
	
    if( search_delay ){
        gate_delay /= 1000; // to ns
    }
    
	/// get input slew
	std::vector<std::string> input_waveform = split( circuit->get_input_waveform(), ' ');
		
	/// initialization for STAGE 1 Simulation
	const static float vdd = this->get_vdd();
	const static float cmax_range = this->get_cmax_range();
	const static float cmin_range = this->get_cmin_range();
	
    float cal_gate_delay, cal_intr_delay;
    
	Waveform i_driver_cat(vdd,ratio);
	Waveform i_driver_1(vdd,ratio);
	Waveform i_driver_2(vdd,ratio);
	
	float tr = std::stof( input_waveform[1] )/1000; 	// to ns
	dctk::PiModel interconnect = circuit->get_load_interconnect();
	float cn = interconnect.get_cnear();		// ff
	float cf = interconnect.get_cfar();			// ff
	float cr1 = receiver_1.get_cap(3,3);		// pf
	float cr2 = receiver_2.get_cap(3,3);		// pf
	
	float ctotal = ( cn + cf + (cr1+cr2)/2*1000 )/1000; 	// to pf
	
    /// out of Lib range search
    float cmax = ctotal*cmax_range;
    float cmin = ctotal*cmin_range;
    
    if( search_delay ){
        cmax = intr_delay*1.5;	// if searching for delay, use extra wide range
        cmin = intr_delay*0.5;	// in searching mode intr_delay = old ceff.
    }
	
	float ceff1 = (cmax+cmin)/2; 						   	// pf
	float ceff2 = (ctotal - ratio*ceff1)/(1.0-ratio) ;		// pf
    float best_ceff = ceff1, err_delay = gate_delay;
	    
	float ref_time = 0.5*tr ;		// asu_exp waveform
	if( input_waveform[0] == "ramp" ) ref_time = 0.5*tr/4*5 - tr/8;
	float n1_time, n2_time;			// time of reaching waveform's delay trip point
	float n1_50_time, n2_50_time;	// time of reaching 50%
	cal_gate_delay = 0;
	cal_intr_delay = 0;
	
	int max_ite = this->get_max_ite_s1();	// end condition if ite reaches max_ite
    if( search_delay ) max_ite = 100;
    
	int ite = 0;

	/// STAGE 1
	for(ite=0; ite<max_ite; ++ite){
        
		if( search_delay ) {
            ceff1 = cmin + (cmax-cmin)*ite/max_ite; //std::sqrt(cmax*cmin);   // pf
            ceff2 = (ctotal - ratio*ceff1)/(1.0-ratio) ;		// pf
        }
    
		/// get 2 waveforms
		float time_step = current.get_waveform(i_driver_1,0,tr,ceff1);	 // get_waveform returns a proper step size for each current waveform.
		current.get_waveform(i_driver_2,time_step,tr,ceff2);
		/// concat 2 waveforms
		i_driver_cat.concat(i_driver_1,i_driver_2);
		
#ifdef DELAY_CAL_DEBUGGING_MODE
		if(verbose>=2) std::cout << "trans = " << tr << std::endl;
		if(verbose>=2) std::cout << "ceff1 = " << ceff1 << std::endl;
		if(verbose>=2) std::cout << "ceff2 = " << ceff2 << std::endl;
		if(verbose>=2) std::cout << "cmax = " << cmax << std::endl;
		if(verbose>=2) std::cout << "cmin = " << cmin << std::endl;
		if(verbose>=2) std::cout << "i_driver_1 cap = " << std::abs(i_driver_1.total_charge()/vdd) << std::endl;
		if(verbose>=2) std::cout << "i_driver_2 cap = " << std::abs(i_driver_2.total_charge()/vdd) << std::endl;
#endif

		circuit_sim( circuit, i_driver_cat, cr1, cr2, n1_time, n2_time, n1_50_time, n2_50_time, trans_time, (verbose>1) );

		cal_gate_delay = n1_50_time - ref_time;
		cal_intr_delay = n2_50_time - n1_50_time;
		
		
#ifdef DELAY_CAL_DEBUGGING_MODE
		if(verbose>=2) std::cout << "ref time = " << ref_time << std::endl;
		if(verbose>=2) std::cout << "n1 50 time = " << n1_50_time << std::endl;
		if(verbose>=2) std::cout << "n2 50 time = " << n2_50_time << std::endl;
		if(verbose>=2) std::cout << "gate delay = " << cal_gate_delay << std::endl;
		if(verbose>=2) std::cout << "intr delay = " << cal_intr_delay << std::endl;
		if(verbose>=2) std::cout << "--------" << std::endl;
#endif
				
		/// reverse lookup for new ceff1
		float new_ceff = current.lookup_td(n1_time, tr, time_step, ratio);
		
#ifdef DELAY_CAL_DEBUGGING_MODE
		if(verbose>=2) std::cout << "new_ceff = " << new_ceff << std::endl;
#endif

		if( search_delay ){
            
            /// given golden delay, search for ceff1
            if( std::abs(cal_gate_delay-gate_delay) < err_delay ){
                err_delay = std::abs(cal_gate_delay-gate_delay);
                best_ceff = ceff1;
            }
			
        }else{
            
            /// normal mode stage 1
            if(new_ceff==ceff1){
                break;			
            }else if(new_ceff<ceff1){
                cmax = ceff1;
            }else{
                cmin = ceff1;
            }
            
            /// calculate new ceff1, ceff2
            ceff1 = (cmax+cmin)/2;
            ceff2 = (ctotal - ratio*ceff1)/(1.0-ratio);
            
		}
		        
	
#ifdef DELAY_CAL_DEBUGGING_MODE
		float err = new_ceff-ceff1;
		/// output waveform to file
		std::string file_name = "i_driver_" + circuit->get_name()+ "_" + ((current.get_polarity()==POS)?"RISE_":"FALL_") + std::to_string(ite) + ".it";
		if(verbose>=2) i_driver_cat.dump_gnu(file_name);
		if(verbose>=2) std::cout << "err = " << err << std::endl;
	} // end of binary search
	
	if(verbose>=1) std::cout << "n1_time = " << n1_time << " n2_time = " << n2_time << std::endl;
	if(verbose>=1) std::cout << "ceff1 = " << ceff1 << " ceff2 = " << ceff2  << std::endl;
	if(verbose>=2) std::cout << std::endl;
#else	
	} // end of binary search

#endif
    
	
    if( search_delay ) return best_ceff;
    
    gate_delay = cal_gate_delay;
    intr_delay = cal_intr_delay;
	return ceff1;
}


void dump_all_waveform(CurrentTable* current, const std::string& file_name, float time_step){
	
		for(int i=0;i<7;++i){
			for(int j=0;j<7;++j){
				
				Waveform i_driver(current->get_vdd(),current->get_delay_trip());
				std::string ofile_name = "dump_waveform/" + file_name + "_" + std::to_string(i) + "_" + std::to_string(j) + ".it";
				
				current->get_waveform(i_driver, i, j);
				i_driver.dump_gnu(ofile_name);
			}
			
		}
	
	
}

float DelayCalculator::compute_delay_s2(
            CurrentTable& current, ReceiverTable& receiver_1, ReceiverTable& receiver_2, dctk::Circuit* circuit,
            float& gate_delay, float& intr_delay, float& trans_time, float ratio, float ceff1_in
){
    /// Stage 2
    // given ceff1, calculate delay/slew 
    
	/// get input slew
	std::vector<std::string> input_waveform = split( circuit->get_input_waveform(), ' ');
		
	/// initialization for STAGE 1 Simulation
	const static float vdd = this->get_vdd();
	    
	Waveform i_driver_cat(vdd,ratio);
	Waveform i_driver_1(vdd,ratio);
	Waveform i_driver_2(vdd,ratio);
	
	float tr = std::stof( input_waveform[1] )/1000; 	// to ns
	dctk::PiModel interconnect = circuit->get_load_interconnect();
	float cn = interconnect.get_cnear();		// ff
	float cf = interconnect.get_cfar();			// ff
	float cr1 = receiver_1.get_cap(3,3);		// pf
	float cr2 = receiver_2.get_cap(3,3);		// pf
	
	float ctotal = ( cn + cf + (cr1+cr2)/2*1000 )/1000; 	// to pf
	
	float ceff1 = ceff1_in;    // pf
	float ceff2 = (ctotal - ratio*ceff1)/(1.0-ratio) ;		// pf
	 
	float ref_time = 0.5*tr ;		// asu_exp waveform
	if( input_waveform[0] == "ramp" ) ref_time = 0.5*tr/4*5 - tr/8;
	float n1_time, n2_time;			// time of reaching waveform's delay trip point
	float n1_50_time, n2_50_time;	// time of reaching 50%
    
	/// get 2 waveforms
	float time_step = current.get_waveform(i_driver_1,0,tr,ceff1);	 // get_waveform returns a proper step size for each current waveform.
	current.get_waveform(i_driver_2,time_step,tr,ceff2);
	/// concat 2 waveforms
	i_driver_cat.concat(i_driver_1,i_driver_2);
    
    /// circuit simulation
	circuit_sim( circuit, i_driver_cat, cr1, cr2, n1_time, n2_time, n1_50_time, n2_50_time, trans_time );

	gate_delay = n1_50_time - ref_time;
	intr_delay = n2_50_time - n1_50_time;
    
    return current.lookup_td(n1_time, tr, time_step, ratio);
    
}

/// compute delay for one circuit
void DelayCalculator::compute_delay(dctk::CellLib* cell_lib, dctk::Circuit* circuit){
   
	/// load circuit
    
	std::vector<std::string> driver_pin = split( circuit->get_driver(), '/' );
	std::string input_pin_name = driver_pin[1];	// [0]=Instance name, [1]=input pin, [2]=output pin
	std::string cell_name = circuit->get_driver_celltype();
       
	std::vector<std::string> load_pin = split( circuit->get_load(), '/' );
	std::string load_pin_name = load_pin[1];	// [0]=Instance name, [1]=input pin, [2]=output pin
	std::string load_cell_name = circuit->get_load_celltype();
	
	
    /// load current data rise
    int this_thread = omp_get_thread_num();
    CurrentTable& current_table = current_tables[ this_thread ];
	current_table.load_current_data( get_current_data_rise(cell_name), POS );
	set_current_table(current_table);
	
    
	/// load receiver from lib
	ReceiverTable& rc1_rise = rc1_rise_m[ this_thread ];
	ReceiverTable& rc2_rise = rc2_rise_m[ this_thread ];
	ReceiverTable& rc1_fall = rc1_fall_m[ this_thread ];
	ReceiverTable& rc2_fall = rc2_fall_m[ this_thread ];
	rc1_rise.load_data( find_c1_rise(cell_lib,load_cell_name,load_pin_name) );
	rc2_rise.load_data( find_c2_rise(cell_lib,load_cell_name,load_pin_name) );
	rc1_fall.load_data( find_c1_fall(cell_lib,load_cell_name,load_pin_name) );
	rc2_fall.load_data( find_c2_fall(cell_lib,load_cell_name,load_pin_name) );

    
	/// Stage 1
	float gate_delay_r, intr_delay_r, gate_trans_r;
	float gate_delay_f, intr_delay_f, gate_trans_f;
       
    float ceff_rise = compute_delay_s1(current_table, rc1_rise, rc2_rise, circuit, gate_delay_r, intr_delay_r, gate_trans_r, get_cat_ratio_rise(), 2); //std::cout << std::endl;
       
    /// load current data fall
	current_table.load_current_data( get_current_data_fall(cell_name), NEG );
	float ceff_fall = compute_delay_s1(current_table, rc1_fall, rc2_fall, circuit, gate_delay_f, intr_delay_f, gate_trans_f, get_cat_ratio_fall(), 2);
    
    if(get_output_feature()==1){
        
        std::vector<std::string> input_waveform = split( circuit->get_input_waveform(), ' ');
        float tr = std::stof( input_waveform[1] );     // ps
        dctk::PiModel interconnect = circuit->get_load_interconnect();
        float cn = interconnect.get_cnear();	// ff
        float res= interconnect.get_res();		// ohm
        
        fp() << gate_delay_r*1000 << " "
                << gate_delay_f*1000 << " "
                << gate_trans_r*1000 << " "
                << gate_trans_f*1000 << " "
                << res << " "
                << cn << " "
                << tr << " "
                << res*cn << " "
                << ceff_rise*1000 << " "
                << ceff_fall*1000 << " "
                << res*ceff_rise << " "
                << res*ceff_fall << " "
                << intr_delay_r*1000 << " "
                << intr_delay_f*1000 << " "
                << cell_name << " "
                << load_cell_name << std::endl;
    }

	/// set values for delays
	circuit->set_ccs_rise_delay( gate_delay_r * (cell_lib->get_scale_to_ps()) );
	circuit->set_ccs_fall_delay( gate_delay_f * (cell_lib->get_scale_to_ps()) );
	circuit->set_ccs_rise_slew(  gate_trans_r * (cell_lib->get_scale_to_ps()) );
	circuit->set_ccs_fall_slew(  gate_trans_f * (cell_lib->get_scale_to_ps()) );
       
   
}


int compute_delays(dctk::CellLib* cell_lib, dctk::CircuitPtrVec* circuitMgr, spef::Spef* spef) {
    
    
	std::cout << "\n------------" << std::endl;
	std::cout << "Computing Delays for " << circuitMgr->size()
				<< " circuits." << std::endl;
     
    int numCircuits = circuitMgr->size();
    
    char config_file[] = "dctk.config";
	DelayCalculator DC(config_file);    
    
    omp_set_num_threads(DC.get_num_thread());
	
    /// output training data
    if(DC.get_output_feature()==1){
        DC.fp().open("data.txt");
    }
    
	/// load RC from spef to circuit
    set_interconnect_RC(circuitMgr,spef);
    	
	
    std::cout << "load current data" << std::endl;
    /// load current data from lib
    DC.load_current_data(cell_lib);
    
    std::cout << "compute delays" << std::endl;
    
    #pragma omp parallel 
    {  
    #pragma omp for nowait
    for (int i=0; i < numCircuits; i++) {
		DC.compute_delay( cell_lib, (*circuitMgr)[i] );
    } // end of calculate all circuits
    }
    
    
    if(DC.get_output_feature()==1){
        DC.fp().close();
    }
    
    /// boost to all circuits
    
    if( DC.get_boost_mode()!=0 ){
        DC.booster( circuitMgr );
    }
    

    return 0;
}

