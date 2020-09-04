
#include "waveform.hpp"

void Waveform::calculate_delay(){
	
	if(waveform_type==CURRENT){
		/// calculate delay time for a current waveform
	
		float vt=0, vt_next=0;
		if(polarity==NEG) vt_next=vdd;
		
		int j, tsize = this->size();
        
		for(j=0; j<tsize-1; ++j){
			vt = vt_next;
			vt_next = vt + (_t_data[j+1]-_t_data[j])*(_i_data[j+1]+_i_data[j])/2.0 /ceff;
			
			if( (vt_next-delay_trip_v)*(vt-delay_trip_v) <= 0 ) break;
		}
		
		tdelay = _t_data[j] + (_t_data[j+1]-_t_data[j]) * (delay_trip_v-vt)/(vt_next-vt);
		td_index = j;
	}
	
	
	good_delay=true;
}


void Waveform::concat(Waveform& wave1, Waveform& wave2){
	
	int ti1 = wave1.get_td_index();
	int ti2 = wave2.get_td_index();
	
	float dt = wave1.time_step();
	int tsize = ti1 - ti2 + wave2.size();	// the size of concatenated waveform
	this->resize(tsize,0.0);
	
	
	for(int i=0; i<tsize; ++i){
		_t_data[i] = i*dt;
	}
	for(int i=0; i<ti1; ++i){
		_i_data[i] = wave1.i().at(i);
	}
	int j=ti2;
	for(int i=ti1; i<tsize; ++i){
		_i_data[i] = wave2.i().at(j++);
	}
	
	set_dt(dt);
	set_voltage(wave1.get_vdd(), wave1.get_delay_trip());
	set_waveform_type(wave1.get_waveform_type());
	set_polarity(wave1.get_polarity());
	set_ceff(wave1.get_ceff());
	
}

void Waveform::dump_gnu(const std::string& fileName)
{
    std::ofstream fout;
    fout.open(fileName);
    
    for (int i = 0; i < this->size(); i++) {
        fout << _t_data[i] << " " << _i_data[i] << std::endl;
    }
    fout.close();

#ifdef DELAY_CAL_DEBUGGING_MODE
    std::cout << "dump waveform to " << fileName << std::endl;
#endif
}

void Waveform::trace_current(float tail_50_point){
	/// trace for 4 points on current wavefomr
	/// start point, max/min point, tail 50% point, end point
	
	if(waveform_type==VOLTAGE) return;
	
	float y_max=0, x_max=-1;
	float y_50 =0, x_50 =-1;
	for(int i=0; i<this->size()-1; ++i){
		
		float x = _t_data.at(i);
		float y = _i_data.at(i);
		
		if(i==0){
			this->x1 = x;
			this->y1 = y;	
			index1=i;
		}
		
		float abs_y = std::abs(y);
		if( abs_y > std::abs(y_max) ){
			x_max = x;
			y_max = y;
			y_50 = y_max * tail_50_point;
			index2=i;
		}
		
		float x_1 = _t_data.at(i+1);
		float y_1 = _i_data.at(i+1);
		if( (y_50-y)*(y_50-y_1) <= 0 ){
			/// y[i] <= y_50 <= y[i+1]
			x_50 = x_1 - (x_1-x) * (y_1-y_50) / (y_1-y);
			index3=i+1;
			break;
		}
		
	}
	
	x2 = x_max;
	y2 = y_max;
	x3 = x_50;
	y3 = y_50;
	
	this->resize(this->size()+1,0);
	
	/// add 50% tail point to waveform
	for(int i=this->size()-1; i>index3; --i){
		_t_data.at(i) = _t_data.at(i-1);
		_i_data.at(i) = _i_data.at(i-1);
	}
	_t_data.at(index3)=x3;
	_i_data.at(index3)=y3;
	
	x4 = _t_data.at(this->size()-1);
	y4 = _i_data.at(this->size()-1);
	index4=this->size()-1;
	
}

float Waveform::total_charge(){
	
	float q = 0;
	
	for(int i=0; i<this->size()-1; ++i){
		q += (_i_data[i]+_i_data[i+1]) * (_t_data[i+1]-_t_data[i]) / 2;
	}
	return q;
}
void Waveform::get_projection(Waveform* w1, Waveform* w2, float ratio){
	
	/// calculated projection target points
	this->x1 = w1->x1 + ratio * (w2->x1 - w1->x1);
	this->y1 = w1->y1 + ratio * (w2->y1 - w1->y1);
	this->x2 = w1->x2 + ratio * (w2->x2 - w1->x2);
	this->y2 = w1->y2 + ratio * (w2->y2 - w1->y2);
	this->x3 = w1->x3 + ratio * (w2->x3 - w1->x3);
	this->y3 = w1->y3 + ratio * (w2->y3 - w1->y3);
	this->x4 = w1->x4 + ratio * (w2->x4 - w1->x4);
	this->y4 = w1->y4 + ratio * (w2->y4 - w1->y4);
	
	this->index1 = w1->index1;
	this->index2 = w1->index2;
	this->index3 = w1->index3;
	this->index4 = w1->index4;
	
	this->resize(w1->size(),0.0);
	
	float a,b,c,d;
	
	/// reconstruct point 1 to point 2
	solve_4by4_matrix(&a,&b,&c,&d,w1->x1,0,w1->x2,w1->y2,x1,0,x2,y2);	// force waveform start at y=0;
	for(int i=0; i<index2; ++i){
		_t_data[i] = a*(w1->_t_data[i]) + b*(w1->_i_data[i]);
		_i_data[i] = c*(w1->_t_data[i]) + d*(w1->_i_data[i]);
	}
	/// reconstruct point 2 to point 3
	solve_4by4_matrix(&a,&b,&c,&d,w1->x2,w1->y2,w1->x3,w1->y3,x2,y2,x3,y3);	// using tail 50 point
	for(int i=index2; i<index3; ++i){
		_t_data[i] = a*(w1->_t_data[i]) + b*(w1->_i_data[i]);
		_i_data[i] = c*(w1->_t_data[i]) + d*(w1->_i_data[i]);
	}
	/// reconstruct point 3 to point 4
	for(int i=index3; i<=index4; ++i){
		_t_data[i] = a*(w1->_t_data[i]) + b*(w1->_i_data[i]);
		_i_data[i] = c*(w1->_t_data[i]) + d*(w1->_i_data[i]);
	}
	
}

void Waveform::solve_4by4_matrix(float* a, float* b, float* c, float* d,
					   float x0, float y0, float x1, float y1,
					   float xt0, float yt0, float xt1, float yt1)
{
	float inv_det = 1.0/(x0*y1 - x1*y0);
	*a = inv_det * (y1*xt0 - y0*xt1);
	*b = inv_det * (-x1*xt0+ x0*xt1);
	*c = inv_det * (y1*yt0 - y0*yt1);
	*d = inv_det * (-x1*yt0+ x0*yt1);
}
