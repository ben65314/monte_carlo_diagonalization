#include "paramReader.h"

typedef StatesR_T<sType,vType> arrType;

int main(int argc, char *argv[]){
	//Args
	if(argc<2||argc>3){std::cout<<"Number of parameters invalid!\nYou should execute TD_solver with this format:\n TD_solver {paramFile} *{verbose}\nWhere * are optionnal arguments\n\nVerbose values:\n* =0 :(default) Minimal prints\n* >0 : Time steps\n* >4 : Lanczos steps\n* >99 : All prints\n"; exit(0);}
	else if (argc == 3){verbose = std::stoi(argv[2]);}

	//Parallelisation
	//omp_set_num_threads(NUM_THREADS_USED);
	//Initialize seed
	srand(clock());
	
	//Start time reading and prep time
	auto step1 = std::chrono::high_resolution_clock::now();

	//Creation of the outFile
	char cwd[PATH_MAX];
	if (getcwd(cwd, sizeof(cwd)) != NULL) {} 
	else {perror("getcwd() error");}
	
	std::string outFileName = "/out.txt";
	const std::string outFileDir = cwd + outFileName;
	

	//Read file of parameters
	justManyVariables jMV = readParameters("/"+(std::string)argv[1]);

	if (verbose > 0) std::cout<<"Bytes for a state:"<<sizeof(sType)<<"\tMax threads:"<<NUM_THREADS_USED<<std::endl;
	//Parallel implementation
	if(jMV.hubP.n_sites > 10){omp_set_num_threads(NUM_THREADS_USED);}
	else {omp_set_num_threads(1);}
	
	//Variable declaring from the read file
	greenParam gP = jMV.gP;
	//gP.print();
	//jMV.sP.print();
	//jMV.hubP.print();

	//Sampling size verification
	sType gotSubSpaceLen = 0;
	//std::cout<<INIT STATE
	Electrons elec = transformNSz(jMV.hubP.nElec,jMV.hubP.totSpin);

	if(jMV.sP.reticle == 0) jMV.sP.reticle = jMV.sP.samplingSize;

	std::string fundStateString;
	//Start computing 
	auto step2 = std::chrono::high_resolution_clock::now();
	arrType MH_Block;
	MH_Block.set_hubbard_parameters(jMV.hubP);
	MH_Block.set_sampling_parameters(jMV.sP);
	MH_Block.electrons = elec;

	//Add the initial states
	for (uLong i = 0; i < jMV.sP.initState.size(); i++) {
		MH_Block.add(jMV.sP.initState.at(i));
	}

	auto step2_1 = std::chrono::high_resolution_clock::now();
	if(verbose > 0){std::cout << "Step 1:Choosing States..."; std::cout.flush();}	

	//Sampling methods
	//MH_Block.sampling_MH();
	MH_Block.sampling_least_energy();
	//MH_Block.showAllStates();

	auto step2_2 = std::chrono::high_resolution_clock::now();
	if(verbose > 0){std::cout << "(Completed) " << timeFormating(step2_1,step2_2) << std::endl;}
	gotSubSpaceLen = MH_Block.getLength();

	//Finds the fundamental energy of the block
	auto step2_3 = std::chrono::high_resolution_clock::now();
	if(verbose > 0){std::cout << "Step 2:Fundamental state...";std::cout.flush();}
	std::vector<vType> fundState = std::vector<vType>(gotSubSpaceLen,0);
	int deg = 1;
	LanczosSolver<vType,decltype(MH_Block)> LS;
	double fundE = LS.fundEnergy(&fundState, &MH_Block, &deg);
	if(verbose > 99) std::cout<<"ENDING FUND VECTOR"<<std::endl;
	if(verbose > 99) printVector(fundState.data(),gotSubSpaceLen);
	//MH_Block.showAllStates();
	if (MH_Block.sys_sP.fund_tc < 1) { 
		if (verbose > 9) std::cout<<"TRUNCATING"<<std::endl;
		if (verbose > 9) std::cout<<"SIZE BEFORE:"<<MH_Block.getLength()<<std::endl;
		writeStateWithDouble(&fundState,&MH_Block,MH_Block.sys_hubP.n_sites,MH_Block.sys_sP.fund_tc);
		if (verbose > 9) std::cout<<"SIZE AFTER:"<<MH_Block.getLength()<<std::endl;
	}

	//printVector(fundState,gotSubSpaceLen);
	auto step2_4 = std::chrono::high_resolution_clock::now();
	if(verbose > 0){std::cout << "(Completed) " << timeFormating(step2_3,step2_4) << std::endl;}
	
	if(jMV.gP.g_compute){
		auto step2_5 = std::chrono::high_resolution_clock::now();
		if(verbose > 0){std::cout << "Step 3:Green functions..."; std::cout.flush();}
		if (MH_Block.sys_hubP.n_sites < 5 ) {
			computeGreen_long(gP.g_AddedSpin, &fundState, fundE, &MH_Block, gP, deg);
		}
		else { 
			computeQMatrixBandLanczos(gP.g_AddedSpin, fundState.data(), fundE, &MH_Block, gP,deg);
		}
		//computeGreen_long(gP.g_AddedSpin, &fundState, fundE, &MH_Block, gP, deg);
		//computeGreen_continued_fraction(gP.g_AddedSpin,fundState.data(),fundE,&MH_Block,gP,deg);

		auto step2_6 = std::chrono::high_resolution_clock::now();
		if(verbose > 0){std::cout << "(Completed) " << timeFormating(step2_5,step2_6)<<std::endl;}
	}
	fundStateString = writeVector(fundState.data(), gotSubSpaceLen);
	//delete[] fundState;
	
	//Start writing
	auto step3 = std::chrono::high_resolution_clock::now();
	
	//WRITING
	std::time_t end_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	std::string writes = (std::string)"\nRun Parameters: --" + std::ctime(&end_time) + (std::string)"#Sites = "+ toStringP(MH_Block.sys_hubP.n_sites,1)+"\tu = "+toStringP(MH_Block.sys_hubP.u,1)+"\tmu = "+ toStringP(MH_Block.sys_hubP.mu,1) + "\n";
	float per = (float)MH_Block.getLength() / (comb(MH_Block.sys_hubP.n_sites,elec.up) * comb(MH_Block.sys_hubP.n_sites,elec.down));
	writes +="Initial State = " + writeVector(MH_Block.sys_sP.initState.data(),MH_Block.sys_sP.initState.size());
	writes +="\nSampling Size  = " + toStringP(MH_Block.getLength(),0) + "(" + toStringP(per,4) +")\tReticle = "+ toStringP(MH_Block.sys_sP.reticle,1) +"\tBeta = " + toStringP(MH_Block.sys_sP.beta_MH,3) + "\n";
	if(verbose == 2) {writes += "States taken = " + MH_Block.showAllStatesString()+"\n";}


	writes += "\nFund energy : " + toStringP(fundE); 
	writes += "\n----------------------------------------------\n";
	
	//Out and writings in files
	std::ofstream outFile;
	outFile.open(outFileDir,std::ios::app);
	
	std::cout<<writes<<std::endl;
	outFile<<writes;
	
	outFile.close();

	//End time
	auto step4 = std::chrono::high_resolution_clock::now();
	if(verbose > 0){
		std::cout<<"Reading time : "<<timeFormating(step1,step2)<< std::endl;
		std::cout<<"Compute time : "<<timeFormating(step2,step3)<< std::endl;
		std::cout<<"Writing time : "<<timeFormating(step3,step4)<< std::endl;
		std::cout<<"Total time : "<<timeFormating(step1,step4)<< std::endl;
	}
	return gotSubSpaceLen;
}
