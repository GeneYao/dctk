#ifndef _CURRENT_TABLE_HPP
#define _CURRENT_TABLE_HPP

#include "dctk.hpp"
#include "CellLib.hpp"
#include "Circuit.hpp"
#include <parser-spef.hpp>
#include "si2dr_liberty.h"

#include "waveform.hpp"
#include <algorithm>

class CurrentTableData{
    
    float* tran;
    float* ceff;
    std::vector<float>** time;
    std::vector<float>** iout;
	int ceff_index_size;
	int tran_index_size;    
    
public:
	void set_ceff_index_size(int i){
		ceff_index_size = i;
	}
	void set_tran_index_size(int i){
		tran_index_size = i;
	}
    void read_table(const std::vector<liberty_value_data*>& table_in);
    
    float* get_tran() const{
        return tran;
    }
    float* get_ceff() const{
        return ceff;
    }
    std::vector<float>** get_time() const{
        return time;
    }
    std::vector<float>** get_iout() const{
        return iout;
    }
    
    
    CurrentTableData(const std::vector<liberty_value_data*>& table_in){
        set_ceff_index_size(7);
        set_tran_index_size(7);
        
        tran = new float [tran_index_size];
        ceff = new float [ceff_index_size];
        time = new std::vector<float> * [tran_index_size];
        for (int i = 0; i < tran_index_size; i++) {
            time[i] = new std::vector<float> [ceff_index_size];
        }
        iout = new std::vector<float> * [tran_index_size];
        for (int i = 0; i < tran_index_size; i++) {
            iout[i] = new std::vector<float> [ceff_index_size];
        }
	
        read_table(table_in);
    }

    ~CurrentTableData(){
        
        delete []tran;
        delete []ceff;
        for (int i = 0; i < tran_index_size; i++) {
            delete [] time[i];
        }
        for (int i = 0; i < tran_index_size; i++) {
            delete [] iout[i];
        }
        delete [] time;
        delete [] iout;
    }
  
};

class CurrentTable {
    
    public:    
        CurrentTable()
            : _input_t(0), _input_c(0), target_table_index(std::pair<int,int>(0,0)),_polarity(POS),
			vdd(0.7), delay_trip(0.5), tr(0.5), tr_max(1), tr_min(0), cr(0.5), cr_max(1), cr_min(0)
        {
            /// default parameter
            set_ceff_index_size(7);
            set_tran_index_size(7);
            set_num_step(100);
            set_max_ite_t1(10);
            set_max_ite_t2(10);
            set_max_ite_cap(10);
            set_max_ite_lookup(10);
            set_project_tail_point(0.5);
            
        }   
        CurrentTable(const CurrentTableData& data, Polarity pol)
            : _input_t(0), _input_c(0), target_table_index(std::pair<int,int>(0,0)),_polarity(pol),
			vdd(0.7), delay_trip(0.5), tr(0.5), tr_max(1), tr_min(0), cr(0.5), cr_max(1), cr_min(0)
        {
            /// default parameter
            set_ceff_index_size(7);
            set_tran_index_size(7);
            set_num_step(100);
            set_max_ite_t1(10);
            set_max_ite_t2(10);
            set_max_ite_cap(10);
            set_max_ite_lookup(10);
            set_project_tail_point(0.5);
            
            load_current_data(data, pol);
        }        
        
        void load_current_data(const CurrentTableData& data, Polarity pol){
            set_polarity(pol);
            tran = data.get_tran();
            ceff = data.get_ceff();
            time = data.get_time();
            iout = data.get_iout();
        }
        
        void reset_waveforms(){
        
            i_tr_c1.resize(num_step, 0);
            i_t2_c1.resize(num_step, 0);
            i_t1_c2.resize(num_step, 0);
            i_t1_c1.resize(num_step, 0);
            i_t2_c2.resize(num_step, 0);
            i_tr_c2.resize(num_step, 0);
            i_composite.resize(num_step, 0);
            i_upper.resize(num_step, 0);
            i_lower.resize(num_step, 0);
            i_find.resize(num_step, 0);
            
        }
        
        void Search(float input_t, float input_c);
        void interpolate_1D(const std::vector<float>& y,const std::vector<float>& t
				   ,std::vector<float>& qy, std::vector<float>& qt
				   ,float dt);

        void get_waveform(std::vector<float>& waveform, std::vector<float>& time
						 ,float dt, float input_t, float input_c);	
		
		float get_waveform(Waveform& waveform
						 ,float dt, float input_t, float input_c);
						 
		void get_waveform(Waveform& waveform ,int i, int j);
		void waveform_projection(Waveform* i_composite);
		
		/// set methods
		void set_num_step(float s){
			num_step = s;
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
		void set_project_tail_point(float d){
			project_tail_point = d;
		}
		void set_vdd(float d){
			vdd = d;
            i_tr_c1.set_voltage(vdd,delay_trip);
            i_t2_c1.set_voltage(vdd,delay_trip);
            i_t1_c2.set_voltage(vdd,delay_trip);
            i_t1_c1.set_voltage(vdd,delay_trip);
            i_t2_c2.set_voltage(vdd,delay_trip);
            i_tr_c2.set_voltage(vdd,delay_trip);
            i_composite.set_voltage(vdd,delay_trip);
            i_upper.set_voltage(vdd,delay_trip);
            i_lower.set_voltage(vdd,delay_trip);
            i_find.set_voltage(vdd,delay_trip);
		}
		void set_delay_trip(float d){
			delay_trip = d;
            i_tr_c1.set_voltage(vdd,delay_trip);
            i_t2_c1.set_voltage(vdd,delay_trip);
            i_t1_c2.set_voltage(vdd,delay_trip);
            i_t1_c1.set_voltage(vdd,delay_trip);
            i_t2_c2.set_voltage(vdd,delay_trip);
            i_tr_c2.set_voltage(vdd,delay_trip);
            i_composite.set_voltage(vdd,delay_trip);
            i_upper.set_voltage(vdd,delay_trip);
            i_lower.set_voltage(vdd,delay_trip);
            i_find.set_voltage(vdd,delay_trip);
		}
		void set_ceff_index_size(int i){
			ceff_index_size = i;
		}
		void set_tran_index_size(int i){
			tran_index_size = i;
		}
		void set_tr(float d){
			tr = d;
		}
		void set_tr_max(float d){
			tr_max = d;
		}
		void set_tr_min(float d){
			tr_min = d;
		}
		void set_cr(float d){
			cr = d;
		}
		void set_cr_max(float d){
			cr_max = d;
		}
		void set_cr_min(float d){
			cr_min = d;
		}
        void set_polarity(Polarity pol){
            _polarity = pol;
        }
		
		/// get methods
		
		float	get_tran_max(){
			return tran[6];
		}
		float	get_tran_min(){
			return tran[0];
		}
		float	get_tran(int index){
			if(index>=0 && index<=6) return tran[index];
			return -1;
		}		
		float	get_ceff_max(){
			return ceff[6];
		}
		float	get_ceff_min(){
			return ceff[0];
		}
		float	get_ceff(int index){
			if(index>=0 && index<=6) return ceff[index];
			return -1;
		}
		float get_vdd(){
			return vdd;
		}
		float get_delay_trip(){
			return delay_trip;
		}
		
		bool is_tran_in_range(float t){
			return t>=get_tran_min() && t<=get_tran_max();
		}
		bool is_ceff_in_range(float c){
			return c>=get_ceff_min() && c<=get_ceff_max();
		}
		float lookup_td(float td, float tr, float dt=0.001, float ratio=0.5);
		Polarity get_polarity(){
			return _polarity;
		}
		
    private:
        float _input_t;
        float _input_c;
        float* tran;
        float* ceff;
        std::vector<float>** time;
        std::vector<float>** iout;
        std::pair<int, int> target_table_index;
		Polarity _polarity;
		float num_step;
		float project_tail_point;
		int max_ite_t1;
		int max_ite_t2;
		int max_ite_cap;
		int max_ite_lookup;
        Waveform i_tr_c1;
        Waveform i_t2_c1;
        Waveform i_t1_c2;
        Waveform i_t1_c1;
        Waveform i_t2_c2;
        Waveform i_tr_c2;
        Waveform i_composite;
        Waveform i_upper;
        Waveform i_lower;
        Waveform i_find;
		int ceff_index_size;
		int tran_index_size;
		float vdd;
		float delay_trip;
		float tr;
		float tr_max;
		float tr_min;
		float cr;
		float cr_max;
		float cr_min;
        
};
#endif

