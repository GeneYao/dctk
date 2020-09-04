#include <iostream>

#include "dctk.hpp"
#include "CellLib.hpp"
#include "Circuit.hpp"
#include <parser-spef.hpp>
#include <string>
#include <vector>
#include <cmath>
#include "PiModel.hpp"
#include "waveform.hpp"
#include "delay_calculator.hpp"
#include <omp.h>
#include <cassert>

const std::vector<liberty_value_data*>& find_current_rise(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name){
	return cell_lib->get_cell(cell_name)->get_output_pin()->find_arc(input_pin_name)->get_output_current_rise_tables() ;
}

const std::vector<liberty_value_data*>& find_current_fall(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name){
	return cell_lib->get_cell(cell_name)->get_output_pin()->find_arc(input_pin_name)->get_output_current_fall_tables() ;
}

const liberty_value_data* find_c1_rise(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name){
	return cell_lib->get_cell(cell_name)->get_output_pin()->find_arc(input_pin_name)->get_receiver_capacitance1_rise_table();
}
const liberty_value_data* find_c2_rise(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name){
	return cell_lib->get_cell(cell_name)->get_output_pin()->find_arc(input_pin_name)->get_receiver_capacitance2_rise_table();
}
const liberty_value_data* find_c1_fall(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name){
	return cell_lib->get_cell(cell_name)->get_output_pin()->find_arc(input_pin_name)->get_receiver_capacitance1_fall_table();
}
const liberty_value_data* find_c2_fall(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name){
	return cell_lib->get_cell(cell_name)->get_output_pin()->find_arc(input_pin_name)->get_receiver_capacitance2_fall_table();
}

const liberty_value_data* find_nldm_rise(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name){
	return cell_lib->get_cell(cell_name)->get_output_pin()->find_arc(input_pin_name)->get_cell_rise_table();
}
const liberty_value_data* find_nldm_fall(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name){
	return cell_lib->get_cell(cell_name)->get_output_pin()->find_arc(input_pin_name)->get_cell_fall_table();
}

std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

void split(std::vector<std::string>& tokens, const std::string& s, char delimiter)
{
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
}

bool compare_circuit_ptr( dctk::Circuit* lhs, dctk::Circuit* rhs ){
	return lhs->get_name() < rhs->get_name();
}
bool compare_spef_ptr( spef::Net* lhs, spef::Net* rhs ){
	return lhs->name < rhs->name;
}

void set_interconnect_RC(dctk::CircuitPtrVec* circuitMgr, spef::Spef* spef){
    
	std::unordered_map< std::string, spef::Net* > net_names;
	
	for( unsigned i=0; i<spef->nets.size(); ++i ){
		net_names[ spef->nets[i].name ] = &(spef->nets[i]);
	}
	
    #pragma omp parallel 
    {  
    #pragma omp for nowait
	for( unsigned i=0; i<circuitMgr->size(); ++i ){
		
		const std::string& str_c = (*circuitMgr)[i]->get_name();
		spef::Net* n = net_names[ str_c ];
		
		if( n==nullptr ){
			std::cerr << "ERROR circuit: " << str_c << " has no spef annotation found" << std::endl;
			assert( n!=nullptr );
		}
		
		/// write RC information to circuit
		unsigned int index = 0;
		float cnear,res,cfar;
		
		for(const auto& c : n->caps){
			auto& [node1, node2, value] = c;
				
			if(index%2) cfar = value;
			else cnear = value;
				
			index++;
		}
		for(const auto& r : n->ress){
			auto& [node1, node2, value] = r;
			res = value;
			index++;
		}
	
		(*circuitMgr)[i]->set_load_interconnect( std::to_string(cnear) + " " + std::to_string(res) + " " + std::to_string(cfar) );
		
	
	}
	}
	
	/*
	std::vector< dctk::Circuit* > circuit_nets( circuitMgr->size(), nullptr );
	std::vector< spef::Net* > spef_nets( spef->nets.size(), nullptr );
	
	/// get nets
	for( unsigned i=0; i<circuit_nets.size(); ++i ){
		circuit_nets[i] = (*circuitMgr)[i];
	}
	for( unsigned i=0; i<spef_nets.size(); ++i ){
		spef_nets[i] = &(spef->nets[i]);
	}
	
	/// sort them by net_name
	std::sort(circuit_nets.begin(), circuit_nets.end(), compare_circuit_ptr);
	std::sort(spef_nets.begin(), spef_nets.end(), compare_spef_ptr);
	
	/// pairing up circuit and spef
	auto ite_c = circuit_nets.begin();
	auto ite_s = spef_nets.begin();
	while( ite_c!=circuit_nets.end() && ite_s!=spef_nets.end() ){
		
		const std::string& str_c = (*ite_c)->get_name();
		const std::string& str_s = (*ite_s)->name;
		
		if( str_c == str_s ){
			/// write RC information to circuit
			unsigned int index = 0;
			float cnear,res,cfar;
			auto& n = **ite_s;
			
			for(const auto& c : n.caps){
				auto& [node1, node2, value] = c;
					
				if(index%2) cfar = value;
				else cnear = value;
					
				index++;
			}
			for(const auto& r : n.ress){
				auto& [node1, node2, value] = r;
				res = value;
				index++;
			}
	
			(*ite_c)->set_load_interconnect( std::to_string(cnear) + " " + std::to_string(res) + " " + std::to_string(cfar) );
			
			ite_c++;
			ite_s++;
			
		}else if( str_c < str_s ){
			std::cerr << "ERROR circuit: " << str_c << " has no spef annotation found" << std::endl;
			assert(str_c >= str_s);
			ite_c++;
		}else{
			ite_s++;
		}
		
	}
	*/
}


int DelayCalculator::circuit_sim(dctk::Circuit* circuit, Waveform& i_driver, float cr, float cr2, float& n1_time, float& n2_time , float& n1_50_time, float& n2_50_time, float& n1_tran, bool verbose){

    /*
		run simulation for circuit:
		[Driver Cnear Res Cfar Receiver]
		
		default unit of input components are (following DCTK)
				Capacitance:  fF
				Resistance : ohm
				Current    :  mA
				Time       :  ns
	*/
	
	/* simulation parameters */
	float dt = i_driver.dt();
	float vdd = i_driver.get_vdd();
	float target_v = vdd;
	float delay_trip_point = i_driver.get_delay_trip_v();
	float delay_50_point = vdd*0.5;
	float delay_10_point = vdd*0.1;
	float delay_90_point = vdd*0.9;
	unsigned i_size = i_driver.size();
	
	const static unsigned step_limit = get_sim_step_limit();
	const static float tolerance_v = get_vdd() * get_sim_tolerance();	// end simulation if v2 reaches its final value within tolerance
	
#ifdef DELAY_CAL_DEBUGGING_MODE
    if(verbose) std::cout << "  >>Circuit Sim" << std::endl;
#endif

	float vv1[step_limit];
	float vv2[step_limit];
	float ttime[step_limit];

	/// get RC information
	dctk::PiModel interconnect = circuit->get_load_interconnect();
    float cn = interconnect.get_cnear()/1000;
    float res= interconnect.get_res()/1000;
    float cf = interconnect.get_cfar()/1000;
	
	/// aid for small res
	float slew_trip = 1;
	if(res<1){			// res<1k
		slew_trip = 0.5*res + 0.5;
		delay_10_point = vdd * (0.5 - 0.4*slew_trip);
		delay_90_point = vdd * (0.5 + 0.4*slew_trip);
	}
	
	
#ifdef DELAY_CAL_DEBUGGING_MODE
    if(verbose) std::cout << "    Cnear=" << cn << "p res=" << res << "k Cfar=" << cf << "p Cr=" << cr << "p" << std::endl;
#endif
	
	
	float v1, v1_next=0;
	float v2, v2_next=0;
	float t_next=0;
	Polarity i_polarity = i_driver.get_polarity();
	
	/// if waveform is Negative, simulation starts from vdd to 0.
	if( i_polarity == NEG ){
		v1_next = v2_next = vdd;
		target_v = 0;
	}
	
	/// measure the time of nodal voltage reached certain value, e.g. 10% 50% 90%
	float d_v1    = 0;
	float d_v1_50 = 0;
	float d_v1_10 = 0;
	float d_v1_90 = 0;
	float d_v2    = 0;
	float d_v2_50 = 0;
    
#ifdef DELAY_CAL_DEBUGGING_MODE
	if(verbose) std::cout << "    driver current size = " << i_size << std::endl;
	if(verbose){
		std::cout << "    driver current polatiry = ";
		if(i_polarity == POS) std::cout << "POS" << std::endl;
		else std::cout << "NEG" << std::endl;
	}
#endif
	
	/// let Res >= 1 ohms
	if(res<0.001) res=0.001;
	float G=1/res;

	/// Trapezoidal method
	float k11 = 2*(cf+cr)+dt*G, k12 = dt*G;
	float k21 = dt*G,		 	 k22 = 2*cn+dt*G;
	float kk = 1/(k11*k22-k12*k21);
	float g11 = 2*cn-dt*G, g22 = 2*(cf+cr)-dt*G;
	float Is_1=i_driver.i().at(0), Is;
	bool sim_v2_reached_50 = false;
	
	/*
	/// Euler method
	float dtG = dt*G;
	float k11 = cf+cr+dtG, k12 = dtG;
	float k21 = dtG,		k22 = cn+dtG;
	float kk = 1/(k11*k22-k12*k21);
	float Is_1=i_driver.i().at(0), Is;
	*/
	
	unsigned i=0;
	while(i<step_limit-1){

#ifdef DELAY_CAL_DEBUGGING_MODE		
#endif	
		vv1[i] = v1_next;
		vv2[i] = v2_next;
		ttime[i] = t_next;
		
		v1 = v1_next;
		v2 = v2_next;
		t_next = t_next + dt;
	
		Is = Is_1;
		if(i<i_size-1){
			Is_1=i_driver.i().at(i+1);
            if(Is_1==0){
                ++i;
                continue;
            }
		}
		
		/// if v2 reached 50%, change to cr2
		if( !sim_v2_reached_50 ){
			if( (i_polarity==POS && v2>=vdd*0.5) || (i_polarity==NEG && v2<=vdd*0.5) ){
				sim_v2_reached_50=true;
				k11 = 2*(cf+cr2)+dt*G;
				g22 = 2*(cf+cr2)-dt*G;
				kk = 1/(k11*k22-k12*k21);
			}
		}
		
		/// Trapezoidal method
		float I1 = g11*v1 + k12*v2 + dt*(Is_1+Is);
		float I2 = k21*v1 + g22*v2;
		v1_next = kk * (k11*I1 + k12*I2);
		v2_next = kk * (k21*I1 + k22*I2);
		
		
		/*
		/// Euler method
		float I1 = cn*v1 + dt*(Is);
		float I2 = (cf+cr)*v2;
		v1_next = kk * (k11*I1 + k12*I2);
		v2_next = kk * (k21*I1 + k22*I2);
		*/
		
		/// check if v2 reached final value within small tolerance
		if( std::abs(v2_next-target_v) < tolerance_v ){
			break;
		}


		++i;
	}
	
	int total_i = i;
#ifdef DELAY_CAL_DEBUGGING_MODE
	if(total_i==step_limit){
		std::cout << "Warning: <CircuitSim> iteration limit reached" << std::endl;
	}
#endif
	
	bool v1_delay_reached = false;
	bool v1_50_reached = false;
	bool v1_10_reached = false;
	bool v1_90_reached = false;
	bool v2_50_reached = false;
	
	/// Delay, slew measurement
	/// preset measurement as final time, prevent from failing in measurement
	d_v1 	= ttime[total_i-1];
	d_v1_50 = ttime[total_i-1];
	d_v1_10 = ttime[total_i-1];
	d_v1_90 = ttime[total_i-1];
	for(int i=0; i<total_i-1; ++i){
		
		// check if v1 reached delay_trip_point
		if( !v1_delay_reached && (vv1[i]-delay_trip_point)*(vv1[i+1]-delay_trip_point)<0 ){
			v1_delay_reached = true;
			d_v1 = ttime[i] + (ttime[i+1] - ttime[i]) * (delay_trip_point-vv1[i]) / (vv1[i+1]-vv1[i]);
		}
		// check if v1 reached 50%
		if( !v1_50_reached && (vv1[i]-delay_50_point)*(vv1[i+1]-delay_50_point)<0 ){
			v1_50_reached = true;
			d_v1_50 = ttime[i] + (ttime[i+1] - ttime[i]) * (delay_50_point-vv1[i]) / (vv1[i+1]-vv1[i]);
		}
		// check if v1 reached lower slew point (10%)
		if( !v1_10_reached && (vv1[i]-delay_10_point)*(vv1[i+1]-delay_10_point)<0 ){
			v1_10_reached = true;
			d_v1_10 = ttime[i] + (ttime[i+1] - ttime[i]) * (delay_10_point-vv1[i]) / (vv1[i+1]-vv1[i]);
		}
		// check if v1 reached upper slew point (90%)
		if( !v1_90_reached && (vv1[i]-delay_90_point)*(vv1[i+1]-delay_90_point)<0 ){
			v1_90_reached = true;
			d_v1_90 = ttime[i] + (ttime[i+1] - ttime[i]) * (delay_90_point-vv1[i]) / (vv1[i+1]-vv1[i]);
		}
        
		// check if v2 reached 50%
		if( !v2_50_reached && (vv2[i]-delay_50_point)*(vv2[i+1]-delay_50_point)<0 ){
			v2_50_reached = true;
			d_v2_50 = ttime[i] + (ttime[i+1] - ttime[i]) * (delay_50_point-vv2[i]) / (vv2[i+1]-vv2[i]);
		}
		
	}
		
	n1_time = d_v1;
	n2_time = d_v2;
	n1_50_time = d_v1_50;
	n2_50_time = d_v2_50;
	
	n1_tran = (d_v1_90-d_v1_10)/slew_trip;
	if(i_polarity==NEG) n1_tran = -n1_tran;

#ifdef DELAY_CAL_DEBUGGING_MODE	
	
    if(verbose) std::cout << "\n    n1 time = " << n1_time << std::endl;
    if(verbose) std::cout << "    n2 time = " << n2_time << std::endl;
    if(verbose) std::cout << "  >>exit Circuit Sim " << std::endl;
#endif
	
	return i;
}
