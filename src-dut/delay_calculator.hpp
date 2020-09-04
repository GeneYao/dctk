#ifndef DELAY_CALCULATOR_HPP_
#define DELAY_CALCULATOR_HPP_

//#define DELAY_CAL_DEBUGGING_MODE

#include "dctk.hpp"
#include <parser-spef.hpp>
#include "current_table.hpp"
#include "receiver_table.hpp"
#include <string>
#include <fstream>
#include <c_api.h>

class DelayCalculator{
		
	float vdd;
	float cmax_range;
	float cmin_range;
	float sim_tolerance;
	float project_tail_point;
	float cat_ratio_rise;
	float cat_ratio_fall;
	int max_ite_s1;
	int max_ite_t1;
	int max_ite_t2;
	int max_ite_cap;
	int max_ite_lookup;
	int project_number_step;
	int sim_step_limit;
    int num_feature;
    int boost_mode;
    int num_thread;
    int output_feature;
    
    std::map<std::string,int> encoding;     // celltype encoding for ML model
	std::ofstream fp_feature;               // output feature for ML traning
    
    std::map< std::string, int > cell_map;
    std::vector< CurrentTableData* > current_data_rise;
    std::vector< CurrentTableData* > current_data_fall;
    
    CurrentTable* current_tables;
    ReceiverTable* rc1_rise_m;
    ReceiverTable* rc2_rise_m;
    ReceiverTable* rc1_fall_m;
    ReceiverTable* rc2_fall_m;
    
public:
    DelayCalculator(char* filename)
        :current_data_rise(25,nullptr), current_data_fall(25,nullptr)
    {
        /// default parameters
		vdd					= 0.7;
		cmax_range			= 1;
		cmin_range			= 0.5;
		sim_tolerance		= 0.03;
		project_tail_point	= 0.5;
		cat_ratio_rise		= 0.5;
		cat_ratio_fall 		= 0.5;
		max_ite_s1			= 1;
		max_ite_t1			= 7;
		max_ite_t2			= 7;
		max_ite_cap			= 7;
		max_ite_lookup		= 7;
		project_number_step	= 90;
		sim_step_limit		= 1000;
        num_feature         = 58;
        boost_mode          = 1;
        num_thread          = 16;
        output_feature      = 0;
        
        
        // celltype encoding for ML model
        encoding["BUFx10_ASAP7_75t_R"] = 0;
        encoding["BUFx12_ASAP7_75t_R"] = 1;
        encoding["BUFx12f_ASAP7_75t_R"] = 2;
        encoding["BUFx24_ASAP7_75t_R"] = 3;
        encoding["BUFx2_ASAP7_75t_R"] = 4;
        encoding["BUFx3_ASAP7_75t_R"] = 5;
        encoding["BUFx4_ASAP7_75t_R"] = 6;
        encoding["BUFx4f_ASAP7_75t_R"] = 7;
        encoding["BUFx5_ASAP7_75t_R"] = 8;
        encoding["BUFx6f_ASAP7_75t_R"] = 9;
        encoding["BUFx8_ASAP7_75t_R"] = 10;
        encoding["HB1xp67_ASAP7_75t_R"] = 11;
        encoding["HB2xp67_ASAP7_75t_R"] = 12;
        encoding["HB3xp67_ASAP7_75t_R"] = 13;
        encoding["HB4xp67_ASAP7_75t_R"] = 14;
        encoding["INVx11_ASAP7_75t_R"] = 15;
        encoding["INVx13_ASAP7_75t_R"] = 16;
        encoding["INVx1_ASAP7_75t_R"] = 17;
        encoding["INVx2_ASAP7_75t_R"] = 18;
        encoding["INVx3_ASAP7_75t_R"] = 19;
        encoding["INVx4_ASAP7_75t_R"] = 20;
        encoding["INVx5_ASAP7_75t_R"] = 21;
        encoding["INVx6_ASAP7_75t_R"] = 22;
        encoding["INVx8_ASAP7_75t_R"] = 23;
        encoding["INVxp33_ASAP7_75t_R"] = 24;
               
        current_tables = nullptr;
        rc1_rise_m = nullptr;
        rc2_rise_m = nullptr;
        rc1_fall_m = nullptr;
        rc2_fall_m = nullptr;
        
        read_config(filename);
        allocate_for_tables();
	}
	~DelayCalculator(){
        
        for(unsigned i=0; i<current_data_rise.size(); ++i){
            delete current_data_rise[i];
        }
        for(unsigned i=0; i<current_data_fall.size(); ++i){
            delete current_data_fall[i];
        }
        

        if(current_tables != nullptr) delete [] current_tables;
        if(rc1_rise_m != nullptr) delete [] rc1_rise_m;
        if(rc2_rise_m != nullptr) delete [] rc2_rise_m;
        if(rc1_fall_m != nullptr) delete [] rc1_fall_m;
        if(rc2_fall_m != nullptr) delete [] rc2_fall_m;
    }
	void read_config(char* filename);

	/// set methods
	void set_vdd(float v){
		vdd = v;
	}
	void set_cmax_range(float c){
		cmax_range = c;
	}
	void set_cmin_range(float c){
		cmin_range = c;
	}
	void set_sim_tolerance(float t){
		sim_tolerance = t;
	}
	void set_project_tail_point(float p){
		project_tail_point = p;
	}
	void set_cat_ratio_rise(float r){
		cat_ratio_rise = r;
	}
	void set_cat_ratio_fall(float r){
		cat_ratio_fall = r;
	}
	void set_max_ite_s1(int i){
		max_ite_s1 = i;
	}
	void set_max_ite_t1(int i){
		max_ite_t1 = i;
	}
	void set_max_ite_t2(int i){
		max_ite_t2 = i;
	}
	void set_max_ite_cap(int i){
		max_ite_cap = i;
	}
	void set_max_ite_lookup(int i){
		max_ite_lookup = i;
	}
	void set_project_number_step(int i){
		project_number_step = i;
	}
	void set_sim_step_limit(int i){
		sim_step_limit = i;
	}
	void set_boost_mode(int i){
        boost_mode = i;
    }
    void set_num_feature(int i){
        if(i>0) num_feature = i;
    }
    void set_num_thread(int i){
        if(i>0) num_thread = i;
    }
    void set_output_feature(int i){
        output_feature = i;
    }
    
	/// get_methods
	float get_vdd(){
		return vdd;
	}
	float get_cmax_range(){
		return cmax_range;
	}
	float get_cmin_range(){
		return cmin_range;
	}
	float get_sim_tolerance(){
		return sim_tolerance;
	}
	float get_project_tail_point(){
		return project_tail_point;
	}
	float get_cat_ratio_rise(){
		return cat_ratio_rise;
	}
	float get_cat_ratio_fall(){
		return cat_ratio_fall;
	}
	int get_max_ite_s1(){
		return max_ite_s1;
	}
	int get_max_ite_t1(){
		return max_ite_t1;
	}
	int get_max_ite_t2(){
		return max_ite_t2;
	}
	int get_max_ite_cap(){
		return max_ite_cap;
	}
	int get_max_ite_lookup(){
		return max_ite_lookup;
	}
	int get_project_number_step(){
		return project_number_step;
	}
	int get_sim_step_limit(){
		return sim_step_limit;
	}
	int get_boost_mode(){
		return boost_mode;
	}
	int get_num_feature(){
        return num_feature;
    }
    int get_num_thread(){
        return num_thread;
    }
    int get_output_feature(){
        return output_feature;
    }
    
    int get_encode(std::string cell_type){
        return encoding[cell_type];
    }
    std::ofstream& fp(){
        return fp_feature;
    }
    
    void get_boost_feature(float* input_feature, dctk::Circuit* circuit );
    void boost_rise_delay(float* input_feature, dctk::Circuit* circuit);
    void boost_fall_delay(float* input_feature, dctk::Circuit* circuit);
    void boost_rise_slew(float* input_feature, dctk::Circuit* circuit);
    void boost_fall_slew(float* input_feature, dctk::Circuit* circuit);
    void booster(dctk::CircuitPtrVec* circuitMgr);
    
    float compute_delay_s1(CurrentTable& current, ReceiverTable& receiver_1, ReceiverTable& receiver_2, dctk::Circuit* circuit,
						  float& gate_delay, float& intr_delay, float& trans_time, float ratio, int verbose=0, bool search_delay=false);
					
    float compute_delay_s2(CurrentTable& current, ReceiverTable& receiver_1, ReceiverTable& receiver_2, dctk::Circuit* circuit,
						  float& gate_delay, float& intr_delay, float& trans_time, float ratio, float ceff1_in);
                          
    float compute_delay_nldm(ReceiverTable& delay_table, ReceiverTable& receiver_1, ReceiverTable& receiver_2, dctk::Circuit* circuit,
						  float& gate_delay, float& intr_delay, float& trans_time, float ratio, int verbose=0);
					
	// ---- in circuit_sim.cpp ----
	int circuit_sim(dctk::Circuit* circuit, Waveform& i_driver, float cr, float cr2, 
					 float& n1_time, float& n2_time , float& n1_50_time, float& n2_50_time,
					 float& n1_tran, bool verbose=false);
	// ---- in circuit_sim.cpp ----
    
    void load_current_data(dctk::CellLib* cell_lib);
    const CurrentTableData& get_current_data_rise(std::string cell_name){
        return *(current_data_rise[ cell_map[cell_name] ]);
    }
    const CurrentTableData& get_current_data_fall(std::string cell_name){
        return *(current_data_fall[ cell_map[cell_name] ]);
    }
    
	void set_current_table(CurrentTable& current){
		current.set_num_step(project_number_step);
		current.set_max_ite_t1(max_ite_t1);
		current.set_max_ite_t2(max_ite_t2);
		current.set_max_ite_cap(max_ite_cap);
		current.set_max_ite_lookup(max_ite_lookup);
		current.set_project_tail_point(project_tail_point);
		current.set_vdd(vdd);
	}
    
    void compute_delay(dctk::CellLib* cell_lib, dctk::Circuit* circuit);
    void allocate_for_tables(){
        if(current_tables == nullptr) current_tables = new CurrentTable[num_thread];
        if(rc1_rise_m == nullptr) rc1_rise_m = new ReceiverTable[num_thread];
        if(rc2_rise_m == nullptr) rc2_rise_m = new ReceiverTable[num_thread];
        if(rc1_fall_m == nullptr) rc1_fall_m = new ReceiverTable[num_thread];
        if(rc2_fall_m == nullptr) rc2_fall_m = new ReceiverTable[num_thread];
    }
};

int compute_delays(dctk::CellLib* cell_lib, dctk::CircuitPtrVec* circuitMgr, spef::Spef *spef);
float score_delay_rise(float* input);
float score_delay_fall(float* input);
float score_slew_rise (float* input);
float score_slew_fall (float* input);
float score_delay_rise_lin(float* input);
float score_delay_fall_lin(float* input);
float score_slew_rise_lin (float* input);
float score_slew_fall_lin (float* input);
//float score_delays(float* input);
//float score_slews (float* input);

#endif
