function initialize_my_path
if ~isdeployed
current_path = pwd;    
addpath([current_path '/my_discritization/']);
addpath([current_path '/my_numerical_method/']);
addpath([current_path '/my_functions/']);
addpath([current_path '/my_tpf_setup/']);
addpath([current_path '/my_output/']);
addpath([current_path '/my_input_field/']);
addpath([current_path '/my_plot/']);
addpath([current_path '/my_struct/']);
addpath(genpath([current_path '/my_material/']));

end

end