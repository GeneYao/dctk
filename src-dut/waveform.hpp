#ifndef _WAVEFORM_HPP
#define _WAVEFORM_HPP

#include "si2dr_liberty.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
 
enum WaveformType {CURRENT, VOLTAGE};
enum Polarity {POS, NEG};

class Waveform{
	
    std::vector<float> _i_data;
    std::vector<float> _t_data;
	bool good_delay;
	int td_index;
	float _time_step;
	float tdelay;
	float vdd;
	float delay_trip;
	float delay_trip_v;
	float ceff;
	WaveformType waveform_type;
	Polarity polarity;
	
	float x1,y1,x2,y2,x3,y3,x4,y4;
	float index1,index2,index3,index4;
	
	void solve_4by4_matrix(float* a, float* b, float* c, float* d,
					   float x0, float y0, float x1, float y1,
					   float xt0, float yt0, float xt1, float yt1);
	
public:

	Waveform()
	:x1(0), y1(0), x2(0), y2(0), x3(0),y3(0), x4(0),y4(0),
	 index1(0), index2(0), index3(0), index4(0)
	{
		good_delay=false;
		td_index=0;
		_time_step=0;
		tdelay=0;
		vdd=0.7;
		delay_trip=0.5;
		delay_trip_v=vdd*delay_trip;
		ceff=1;
		waveform_type=CURRENT;
		polarity=POS;
	}
	
	// set vdd, delay trip, WaveformType, Polarity
	Waveform(float vd, float trip, WaveformType wt=CURRENT, Polarity po=POS)
		:good_delay(false), td_index(0), _time_step(0), tdelay(0), vdd(vd), delay_trip(trip), ceff(1),
		 waveform_type(wt), polarity(po), x1(0), y1(0), x2(0), y2(0), x3(0),y3(0), x4(0),
			y4(0), index1(0), index2(0), index3(0), index4(0)
	{
		delay_trip_v=vdd*delay_trip;
	}
	
	// set vdd, delay trip, WaveformType, Polarity and initial size
	Waveform(float vd, float trip, WaveformType wt, Polarity po, unsigned size)
		:_i_data(size,0.0), _t_data(size,0.0),
		 good_delay(false), td_index(0), _time_step(0), tdelay(0), vdd(vd), delay_trip(trip), ceff(1),
		 waveform_type(wt), polarity(po), x1(0), y1(0), x2(0), y2(0), x3(0),y3(0), x4(0),
			y4(0), index1(0), index2(0), index3(0), index4(0)
	{
		delay_trip_v=vdd*delay_trip;
	}
	
	void resize(unsigned s, float v){
		_i_data.resize(s,v);
		_t_data.resize(s,v);	
		good_delay=false;
	}
    
    void clear(){
        _i_data.clear();
        _t_data.clear();
		good_delay=false;
    }
    
	// access i,t vector
	// if you intend to set i,t data
	// please use i_data() so that tdelay will be updated later,
	// otherwise use i() to skip delay time calculation.
	std::vector<float>& i_data(){
		good_delay=false;
		return _i_data;
	}
	const std::vector<float>& i(){
		return _i_data;
	}
	std::vector<float>& t_data(){
		good_delay=false;
		return _t_data;
	}
	const std::vector<float>& t(){
		return _t_data;
	}
	
	// set methods
	void set_dt(float dt){
		good_delay=false;
		_time_step=dt;
	}
	void set_tdelay(float td){
		good_delay=true;
		tdelay=td;
	}
	void set_voltage(float vd, float trip){
		good_delay=false;
		vdd=vd;
		delay_trip=trip;
		delay_trip_v=vd*trip;
	}
	void set_waveform_type(WaveformType wt){
		good_delay=false;
		waveform_type=wt;
	}
	void set_polarity(Polarity po){
		good_delay=false;
		polarity=po;
	}
	void set_ceff(float c){
		good_delay=false;
		ceff=c;
	}
    
	// get methods
	int get_td_index(){
		if(!good_delay) calculate_delay();
		return td_index;
	}
	float get_tdelay(){
		if(!good_delay) calculate_delay();
		return tdelay;
	}
	float get_vdd() const{
		return vdd;
	}
	float get_delay_trip() const{
		return delay_trip;
	}
	float get_delay_trip_v() const{
		return delay_trip_v;
	}
	float get_ceff() const{
		return ceff;
	}
	WaveformType get_waveform_type() const{
		return waveform_type;
	}
	Polarity get_polarity() const{
		return polarity;
	}
	int size() const{
		return _t_data.size();		
	}
	float time_step() const{
		return _time_step;
	}
	float dt() const{
		return _time_step;
	}
    
	// member functions
	void calculate_delay();
	void concat(Waveform& wave1, Waveform& wave2);
	void dump_gnu(const std::string& fileName);
	void trace_current(float tail_50_point=0.5);
	void get_projection(Waveform* w1, Waveform* w2, float ratio);
	float total_charge();
};

#endif