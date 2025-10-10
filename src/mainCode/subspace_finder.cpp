#include "paramReader.h"

//Type to use in code
typedef double vType;
typedef StatesR_H<sType,vType> arrType;

int main(int argc, char *argv[]){
	//Args
	if(argc<2||argc>3){
		std::cout<<"Number of parameters invalid!\n"<<
		"You should execute s_breadthSampling with this format:\n"<<
		"s_breadthSampling {paramFile} *{verbose}\n"<<
		"Where * are optionnal arguments\n\nVerbose values:\n"<<
		"* =0 :(default) Minimal prints\n* >0 : Time steps\n"<<
		"* =2 : Prints Array of State\n"<<
		"* =3 Print Means and deviation\n"<<
		"* >5 : All prints\n"; exit(0);}
	else if (argc == 3){verbose = std::stoi(argv[2]);}

	//Initialize seed
	srand(clock());
	
	//Start time reading and prep time
	auto step1 = std::chrono::high_resolution_clock::now();

	//Read file of parameters
	justManyVariables jMV = readParameters("/"+(std::string)argv[1]);
    
    //Present to the user the number of bytes allowed for states and threads
	if (verbose > 0) std::cout << "Bytes for a state:" << sizeof(sType) <<
		"\tMax threads:" << NUM_THREADS_USED << std::endl;

	//Parallel implementation
	if(jMV.hubP.n_sites > 10){
        omp_set_num_threads(NUM_THREADS_USED);
    }
	else {
        omp_set_num_threads(1);
    }

	//Start computing 
	auto step2 = std::chrono::high_resolution_clock::now();
	
    //Minimum energy
	double minEnergy = 1e10;
	Electrons min_block_electrons;
    
    //Sampling parameters
    samplingParam sP;
    sP.beta_MH = 0;
    sP.nHapply = 0; 
    sP.beta_Happly = 0;

	//Test for a number of electrons up
	for (int i = 0; i <= jMV.hubP.n_sites; i++) {
		//Test for a number of elctrons down
		for (int j = i; j <= jMV.hubP.n_sites; j++) {
			//Skips no electron case
            if (i == j && i == 0) continue;

            //Block size
			sP.reticle = sP.samplingSize = comb(jMV.hubP.n_sites, i)
				*comb(jMV.hubP.n_sites,j);

			Electrons elec;
			elec.up = i; elec.down = j;

            //Setting block
			std::vector<sType> state_init = {
				createAntiFerro(jMV.hubP.n_sites, i, j)};
			arrType states_block = arrType(&state_init);
			states_block.set_hubbard_parameters(jMV.hubP);
			states_block.set_sampling_parameters(sP);
			states_block.electrons = elec;

			//Sampling methods
			states_block.sampling_MH();

            //Find ground state
			std::vector<vType> fundState = std::vector<vType>(
				sP.samplingSize,0);
			LanczosSolver<vType,arrType> LS;
			int deg = 1;
			double fundE = LS.fundEnergy(&fundState, 
				                         &states_block, &deg);
            
            //Look for minimum energy
			if (fundE < minEnergy) {
				minEnergy = fundE;
				min_block_electrons = elec;
			}
		}
	}
	auto step3 = std::chrono::high_resolution_clock::now();

    //Information on min_block
	int N = min_block_electrons.up + min_block_electrons.down;
	int S = min_block_electrons.up - min_block_electrons.down;
	std::cout<<"The fundamental block is N"<<N<<"S"<<S<<
		"\nThe fundamental energy is "<<minEnergy<<std::endl;


	//End time
	if (verbose > 0) {
		std::cout << "Reading time : " << timeFormating(step1, step2) 
			<< std::endl;
		std::cout << "Compute time : " << timeFormating(step2, step3) 
			<< std::endl;
		std::cout << "Total time : " << timeFormating(step1, step3) 
			<< std::endl;
	}
	return 0;
}
