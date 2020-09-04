#ifndef _CIRCUIT_SIM_HPP
#define _CIRCUIT_SIM_HPP

void set_interconnect_RC(dctk::CircuitPtrVec* circuitMgr, spef::Spef* spef);

std::vector<std::string> split(const std::string& s, char delimiter);
void split(std::vector<std::string>& tokens, const std::string& s, char delimiter);

const std::vector<liberty_value_data*>& find_current_rise(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name);
const std::vector<liberty_value_data*>& find_current_fall(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name);

const liberty_value_data*   find_c1_rise(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name);
const liberty_value_data*   find_c2_rise(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name);
const liberty_value_data*   find_c1_fall(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name);
const liberty_value_data*   find_c2_fall(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name);
const liberty_value_data* find_nldm_rise(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name);
const liberty_value_data* find_nldm_fall(dctk::CellLib* cell_lib, const std::string& cell_name, const std::string& input_pin_name);

#endif