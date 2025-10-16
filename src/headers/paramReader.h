#pragma once
#include "greenFunctions.h"

struct justManyVariables readParameters(const std::string file) {
	char cwd[PATH_MAX];
	if (getcwd(cwd, sizeof(cwd)) != NULL) {
		//printf("Current working dir: %s\n", cwd);
	} else {
		perror("getcwd() error");
   	}	
	//File Name
	std::string data_file = cwd+file;
	//Variables read
	int N = 0;
	int Sz = 0;
	Electrons elec;
	sType init;
	std::vector<int> sites_pos;
	unsigned char current_site = 0;
	std::vector<int> allowed_jump;
	std::vector<double> jump_energy;
	std::vector<double> mu_array;
	jump_energy.reserve(5);

	float percentage_of_states = 1;
	
	//Struct countaining all the variables that will be passed to the program
	justManyVariables send_info;

	bool use_specified=false;
	std::vector<sType> init_state;
	
	//Recherche du fichier ayant les donnees necessaire
	std::string lines;
	std::fstream my_file;
	my_file.open(data_file,std::ios::in);
	int number_of_allowed_jump = 0;
	if (my_file.is_open())
	{
		while (getline(my_file, lines))
		{
			//Contains the information between brackets
			int start_info = lines.find('{');
			int end_info = lines.find('}');
			//For vector informations
			int startP = lines.find('(');
			int endP = lines.find(')');
			int del1 = lines.find(',');
			int del2 = lines.find(',', del1 + 1);

            //If the current line has brakets {} to read into
			bool has_brackets = (start_info != -1) 
                            && (end_info != -1) 
                            && (end_info - start_info > 1);
            //If the current line has parentheses () for vector to read into
			bool has_vector = (startP != -1) 
                          && (endP != -1) 
                          && (del1 != -1) 
                          && (del1 != -1);

            //Values read
			uLong value_int = 0;
			float value_float = 0;
			std::string value_str = " ";
			std::vector<int> value_vector;

			if(has_brackets) {
				if(has_vector) {
                    //Read vector
					value_vector.push_back(std::stoi(lines.substr(
                                                startP+1, del1-(startP + 1))));
					value_vector.push_back(std::stoi(lines.substr(
                                                del1+1, del2-(del1 + 1))));
					value_vector.push_back(std::stoi(lines.substr(
                                                del2+1, endP-(del2 + 1))));
				}
				else {
                    //Read value int or float
					value_str = lines.substr(start_info + 1, 
                                             end_info - start_info - 1);
					try {value_int = std::stoll(value_str);} catch(...){}
					try {value_float = std::stof(value_str);} catch(...){}
				}
			}

			//Empty line
			if (lines.size()==0) {
			}
			//Comments
			else if (lines.starts_with("#")) {
			}
			else if (lines.find("number_of_sites") != std::string::npos) {
				send_info.hubP.n_sites = value_int;
			}
			else if (lines.find("sampling_space_size") != std::string::npos) {
				send_info.sP.sampling_size = value_int;
            }
			else if (lines.find("U") != std::string::npos) {
				send_info.hubP.u = value_float;
			}
			else if (lines.find("mu") != std::string::npos) {
				send_info.hubP.mu = value_float;
			}
			else if (lines.find("N") != std::string::npos) {
				N = value_int;
			}
			else if (lines.find("Sz") != std::string::npos) {
				Sz = value_int;
			}
		    else if (lines.find("beta") != std::string::npos) {
				send_info.sP.beta_MH = value_float;
			}
			else if (lines.find("Beta_Happly") != std::string::npos) {
				send_info.sP.beta_Happly = value_float;
			}
			else if (lines.find("nHapply") != std::string::npos) {
				send_info.sP.nHapply = value_int;
			}
			else if (lines.find("reticle") != std::string::npos) {
				send_info.sP.reticle = value_int;
			}	
			else if (lines.find("truncated_cutoff") != std::string::npos) {
				send_info.sP.fund_tc = value_float;
			}
			//Location of sites//////////////////////////////////////
			else if (lines.find("site_location") != std::string::npos) {
				sites_pos.insert(sites_pos.end(), value_vector.begin(),
                                value_vector.end());
				current_site++;

				//Adds variable mu for each sites
				int startB = lines.find('{', endP);
				int endB = lines.find('}', startB);
				if (startB > -1 || endB > -1) {
					double tempo = std::stof(lines.substr(startB + 1, 
                                                    endB - (startB + 1)));
					mu_array.push_back(tempo);
				}
			}
			//Allowed Jumps//////////////////////////////////////////////
			else if (lines.find("allow_jump") != std::string::npos) {
				allowed_jump.insert(allowed_jump.end(), value_vector.begin(),
                                   value_vector.end());
				number_of_allowed_jump++;
        
                //Adds the energy according to the jump
				int startB = lines.find('{', endP);
				int endB = lines.find('}', startB);
				double tempo = std::stof(lines.substr(startB + 1,
                                         endB - (startB + 1)));
                jump_energy.push_back(tempo);
			}
			//Get initial s states
			else if (lines.find("initial_states") != std::string::npos){
				std::vector<int> coma_index;

				//Looks for all the commas and brakets
				coma_index.push_back(start_info);
				int i = start_info;
				while(i>0){
					i = lines.find(',',i+1);
					if(i<0){
						break;
					}
					coma_index.push_back(i);
					
				}
				coma_index.push_back(end_info);

				//What if there are no commas
				if(coma_index.size()>2 || ((end_info - start_info)>1)){
					use_specified =true;
				}
				else{
					continue;
				}
				
				//Extract states
				std::vector<sType> temp_initial_state_num;
				for(unsigned int i = 0; i < coma_index.size()-1; i++){
					std::string a = lines.substr(coma_index.at(i) + 1,
                                  coma_index.at(i+1) - (coma_index.at(i) + 1));
					temp_initial_state_num.push_back(std::strtoll(a.c_str(),
                                                               NULL, 0));

				}	
				init_state = temp_initial_state_num;
			}
			else if (lines.find("added_site") != std::string::npos) {
				send_info.gP.g_added_site = value_int;
			}
			else if (lines.find("added_spin") != std::string::npos) {
				send_info.gP.g_added_spin = value_int;
			}
			else if (lines.find("eta") != std::string::npos) {
				send_info.gP.g_eta = value_float;
			}
			else if (lines.find("add_all") != std::string::npos) {
				send_info.gP.g_add_all = (value_str == "YES");
			}
			else if (lines.find("average") != std::string::npos) {
				send_info.gP.g_average = (value_str == "YES");
			}
			else if (lines.find("compute_green") != std::string::npos) {
				send_info.gP.g_compute = (value_str == "YES");
			}
		}
		my_file.close();
	}
	else { 
        std::cout << "Unable to open the parameter file named : "
            << file << std::endl; 
        exit(1);
    }
	
	//Verification of entered parameters
	//n_sites = 0
	if(send_info.hubP.n_sites == 0) {
        std::cout << "ERROR : The number of sites can't be zero" << std::endl;
        exit(1);
    }
	// N > 2*sites
	if (N > send_info.hubP.n_sites * 2) {
        std::cout << "ERROR : Too many electrons\n" 
            << "You have put more electrons than twice the number of sites" 
            << std::endl; 
        exit(1);
    }
	// sites != site position entered
	if(send_info.hubP.n_sites != current_site){
		std::cout << "ERROR : The number of sites entered doesn't match the" 
            << " number of sites location" << std::endl;
		exit(1);
	}
	//Percentage not between 0 and 1
	if (percentage_of_states < 0 || percentage_of_states > 1) {
		std::cout << "ERROR : Percentage of states R entered" 
            << " not between 0 and 1" << std::endl;
		exit(1);
	}
	
	float up_f = (N+Sz) / 2;
	float down_f = (N-Sz) / 2;
	if ((up_f != (int)up_f) || (down_f != (int)down_f)) {
		std::cout << "N and Sz are incompatible!" << std::endl;	
		exit(1);
	}
	elec.up = up_f;
	elec.down = down_f;
	
	//Creating tMatrix
	double* t_matrix = new double[send_info.hubP.n_sites 
                                            * send_info.hubP.n_sites]();
	for (unsigned char i = 0; i < send_info.hubP.n_sites; i++) {	
		for (unsigned char j = i; j < send_info.hubP.n_sites; j++) {
			if (i == j) {
				t_matrix[i*send_info.hubP.n_sites+j]=0;
				continue;
			}
			int jumpX = sites_pos.at(i*3 + 0) - sites_pos.at(j*3 + 0); 	
			int jumpY = sites_pos.at(i*3 + 1) - sites_pos.at(j*3 + 1); 
			int jumpZ = sites_pos.at(i*3 + 2) - sites_pos.at(j*3 + 2);

			for (int k = 0; k < number_of_allowed_jump; k++) {
				bool good_jump =(jumpX == allowed_jump.at(k*3 + 0) 
                              && jumpY == allowed_jump.at(k*3 + 1) 
                              && jumpZ == allowed_jump.at(k*3 + 2))
							 || (jumpX == -allowed_jump.at(k*3 + 0) 
                              && jumpY == -allowed_jump.at(k*3 + 1) 
                              && jumpZ == -allowed_jump.at(k*3 + 2));
				if (good_jump) {
					t_matrix[i*send_info.hubP.n_sites + j] = jump_energy.at(k);
					t_matrix[j*send_info.hubP.n_sites + i] = jump_energy.at(k);
				}
			}
		}
	}
	
	send_info.hubP.t_matrix = std::vector<double>(
        t_matrix, t_matrix + (send_info.hubP.n_sites*send_info.hubP.n_sites));
	delete[] t_matrix;

	init = create_anti_ferro(send_info.hubP.n_sites, up_f, down_f);

	std::vector<sType> temp;
	temp.push_back(init);
	
	//Using entered or calculated states
	if (percentage_of_states > 1) percentage_of_states = 1;
	if (use_specified) send_info.sP.init_state = init_state;
	else send_info.sP.init_state = temp;
	
	
	send_info.hubP.N_e = N;
	send_info.hubP.S_z = Sz;

	if(send_info.sP.sampling_size <=0) send_info.sP.sampling_size = 1;

	return send_info;
	
}
