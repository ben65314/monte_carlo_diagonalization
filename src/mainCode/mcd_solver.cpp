#include "paramReader.h"

typedef StatesR_T<sType,vType> arrType;

int main(int argc, char *argv[]){
	//Args
	if (argc<2 || argc>3) {
        std::cout << "Number of parameters invalid!\nYou should execute \n" 
                  << "TD_solver with this format:\n TD_solver {paramFile} "
                  << "*{verbose}\nWhere * are optionnal arguments\n\nVerbose "
                  << "values:\n* =0 :(default) Minimal prints\n* >0 : Time " 
                  << "steps\n* >4 : Lanczos steps\n* >99 : All prints\n"; 
        exit(0);
    }
	else if (argc == 3) verbose = std::stoi(argv[2]);

	//Initialize seed
	srand(clock());
	
	//Start time reading and prep time
	auto step1 = std::chrono::high_resolution_clock::now();

	//Creation of the outFile
	char cwd[PATH_MAX];
	if (getcwd(cwd, sizeof(cwd)) != NULL) {}
	else perror("getcwd() error");
	
	std::string out_file_name = "/out.txt";
	const std::string out_file_dir = cwd + out_file_name;

	//Read file of parameters
	justManyVariables jMV = readParameters("/" + (std::string)argv[1]);

	if (verbose > 0) {
        std::cout << "Bytes for a state:" << sizeof(sType) << "\tMax threads:"
                  << NUM_THREADS_USED << std::endl;
    }
	//Parallel implementation
    if (jMV.hubP.n_sites > 10) omp_set_num_threads(NUM_THREADS_USED);
	else omp_set_num_threads(1);
	
	//Variable declaring from the read file
	greenParam gP = jMV.gP;

	//Sampling size verification
	sType got_subspace_len = 0;
	Electrons elec = transform_NSz(jMV.hubP.N_e, jMV.hubP.S_z);

	if (jMV.sP.reticle == 0) jMV.sP.reticle = jMV.sP.sampling_size;

	std::string fund_state_string;
	//Start computing 
	auto step2 = std::chrono::high_resolution_clock::now();
	arrType MH_Block;
	MH_Block.set_hubbard_parameters(jMV.hubP);
	MH_Block.set_sampling_parameters(jMV.sP);
	MH_Block.electrons = elec;

	//Add the initial states
	for (uLong i = 0; i < jMV.sP.init_state.size(); i++) {
		MH_Block.add(jMV.sP.init_state.at(i));
	}

	auto step2_1 = std::chrono::high_resolution_clock::now();
	if (verbose > 0) std::cout << "Step 1:Choosing States...";	

	//Sampling methods
	//MH_Block.sampling_MH();
	MH_Block.sampling_least_energy();

	auto step2_2 = std::chrono::high_resolution_clock::now();
	if (verbose > 0) {
        std::cout << "(Completed)" << time_formating(step2_1, step2_2) << '\n';
    }
	got_subspace_len = MH_Block.get_length();

	//Finds the fundamental energy of the block
	auto step2_3 = std::chrono::high_resolution_clock::now();
	if (verbose > 0) std::cout << "Step 2:Fundamental state...";

	std::vector<vType> fund_state = std::vector<vType>(got_subspace_len, 0);
	int deg = 1;
	LanczosSolver<vType,decltype(MH_Block)> LS;
	double fundE = LS.fund_energy(&fund_state, &MH_Block, &deg);

	if (verbose > 99) std::cout << "ENDING FUND VECTOR" << std::endl;
	if (verbose > 99) print_vector(fund_state.data(), got_subspace_len);
	//MH_Block.show_all_states();
	if (MH_Block.sys_sP.fund_tc < 1) { 
		if (verbose > 9) {
            std::cout << "TRUNCATING" << std::endl;
		    std::cout << "SIZE BEFORE:" << MH_Block.get_length() << std::endl;
        }
		write_state_with_double(&fund_state, &MH_Block,
                                MH_Block.sys_hubP.n_sites,
                                MH_Block.sys_sP.fund_tc);
		if (verbose > 9) {
            std::cout << "SIZE AFTER: " << MH_Block.get_length() << std::endl;
        }
	}

	auto step2_4 = std::chrono::high_resolution_clock::now();
	if (verbose > 0) {
        std::cout << "(Completed)" << time_formating(step2_3, step2_4) << '\n';
    }
	
	if (jMV.gP.g_compute) {
		auto step2_5 = std::chrono::high_resolution_clock::now();
		if (verbose > 0) std::cout << "Step 3:Green functions...";

		if (MH_Block.sys_hubP.n_sites < 5 ) {
			compute_green_long(gP.g_added_spin, &fund_state, fundE, 
                                &MH_Block, gP, deg);
		}
		else { 
			compute_q_matrix_band_lanczos(gP.g_added_spin, fund_state.data(), 
                                            fundE, &MH_Block, gP, deg);
		}

		auto step2_6 = std::chrono::high_resolution_clock::now();
		if (verbose > 0) {
            std::cout << "(Completed) " << time_formating(step2_5, step2_6) 
                      << std::endl;
        }
	}
	fund_state_string = write_vector(fund_state.data(), got_subspace_len);
	//delete[] fundState;
	
	//Start writing
	auto step3 = std::chrono::high_resolution_clock::now();
	
	//WRITING
	std::time_t end_time = std::chrono::system_clock::to_time_t(
                                        std::chrono::system_clock::now());

	std::string writes = (std::string)"\nRun Parameters: --" 
                        + std::ctime(&end_time) + (std::string)"#Sites = "
                        + to_string_p(MH_Block.sys_hubP.n_sites, 1) + "\tu = "
                        + to_string_p(MH_Block.sys_hubP.u, 1) + "\tmu = "
                        + to_string_p(MH_Block.sys_hubP.mu, 1) + "\n";

	float per = (float)MH_Block.get_length() 
                / (comb(MH_Block.sys_hubP.n_sites,elec.up) 
                    * comb(MH_Block.sys_hubP.n_sites,elec.down));

	writes += "Initial State = " 
            + write_vector(MH_Block.sys_sP.init_state.data(), 
                           MH_Block.sys_sP.init_state.size());

	writes += "\nSampling Size  = " 
            + to_string_p(MH_Block.get_length(), 0) + "(" + to_string_p(per, 4) 
            + ")\tReticle = " + to_string_p(MH_Block.sys_sP.reticle, 1) 
            + "\tBeta = " + to_string_p(MH_Block.sys_sP.beta_MH, 3) + "\n";

	if (verbose == 2) {
        writes += "States taken = " + MH_Block.show_all_states_string() + "\n";
    }


	writes += "\nFund energy : " + to_string_p(fundE); 
	writes += "\n----------------------------------------------\n";
	
	//Out and writings in files
	std::ofstream out_file;
	out_file.open(out_file_dir, std::ios::app);
	
	std::cout << writes << std::endl;
	out_file << writes;
	
	out_file.close();

	//End time
	auto step4 = std::chrono::high_resolution_clock::now();
	if(verbose > 0){
		std::cout << "Reading time : " << time_formating(step1, step2) << '\n'
		          << "Compute time : " << time_formating(step2, step3) << '\n'
		          << "Writing time : " << time_formating(step3, step4) << '\n'
		          << "Total time : " << time_formating(step1, step4) << '\n';
	}
	return got_subspace_len;
}
