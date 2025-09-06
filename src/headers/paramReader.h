#pragma once
#include "greenFunctions.h"

//The file name can be change here if you have multiple preset
//const std::string fileRead= "/paramFileCpp.txt";

struct justManyVariables readParameters(const std::string fileRead) {
	char cwd[PATH_MAX];
	if (getcwd(cwd, sizeof(cwd)) != NULL) {
		
		//printf("Current working dir: %s\n", cwd);
	} else {
		perror("getcwd() error");
   	}	
	//File Name
	std::string data_file = cwd+fileRead;
	//Variables read
	int N = 0;
	int Sz = 0;
	Electrons elec;
	sType init;
	std::vector<int> sitesPos;
	unsigned char currentSite = 0;
	std::vector<int> allowedJump;
	std::vector<double> jumpEnergy;
	std::vector<double> mu_array;
	jumpEnergy.reserve(5);

	float percentageOfStates = 1;
	
	//Struct countaining all the variables that will be passed to the program
	justManyVariables sendInfo;

	bool useSpecifiedR=false;
	std::vector<sType> initState;
	
	//Recherche du fichier ayant les donnees necessaire
	std::string lines;
	std::fstream myFile;
	myFile.open(data_file,std::ios::in);
	int numberOfAllowedJump = 0;
	if (myFile.is_open())
	{
		while (getline(myFile, lines))
		{
			//Contains the information between brackets
			int startInfo = lines.find('{');
			int endInfo = lines.find('}');
			//For vector informations
			int startP = lines.find('(');
			int endP = lines.find(')');
			int del1 = lines.find(',');
			int del2 = lines.find(',', del1 + 1);

			bool hasBrackets = (startInfo != -1) && (endInfo != -1) && (endInfo - startInfo > 1);
			bool hasVector = (startP != -1) && (endP != -1) && (del1 != -1) && (del1 != -1);

			uLong value_int = 0;
			float value_float = 0;
			std::string value_str = " ";
			std::vector<int> value_vector;

			if(hasBrackets) {
				if(hasVector) {
					value_vector.push_back(std::stoi(lines.substr(startP+1, del1-(startP + 1))));
					value_vector.push_back(std::stoi(lines.substr(del1+1,del2-(del1 + 1))));
					value_vector.push_back(std::stoi(lines.substr(del2+1, endP-(del2 + 1))));
				}
				else {
					value_str = lines.substr(startInfo + 1, endInfo - startInfo - 1);
					try {value_int = std::stoll(value_str);} catch(...){}
					try {value_float = std::stof(value_str);} catch(...){}
				}
			}

			//std::cout<<lines<<std::endl;
			//Empty line
			if (lines.size()==0){
			}
			//Comments
			else if (lines.starts_with("#")) {
			}
			else if (lines.find("number_of_sites") != std::string::npos) {
				sendInfo.hubP.n_sites = value_int;
			}
			else if (lines.find("sampling_space_size") != std::string::npos) {
				sendInfo.sP.samplingSize = value_int;
            }
			else if (lines.find("U") != std::string::npos) {
				sendInfo.hubP.u = value_float;
			}
			else if (lines.find("mu") != std::string::npos) {
				sendInfo.hubP.mu = value_float;
			}
			else if (lines.find("N") != std::string::npos) {
				N = value_int;
			}
			else if (lines.find("Sz") != std::string::npos) {
				Sz = value_int;
			}
		    else if (lines.find("beta") != std::string::npos) {
				sendInfo.sP.beta_MH = value_float;
			}
			else if (lines.find("Beta_Happly") != std::string::npos) {
				sendInfo.sP.beta_Happly = value_float;
			}
			else if (lines.find("nHapply") != std::string::npos) {
				sendInfo.sP.nHapply = value_int;
			}
			else if (lines.find("reticle") != std::string::npos) {
				sendInfo.sP.reticle = value_int;
			}	
			else if (lines.find("truncated_cutoff") != std::string::npos) {
				sendInfo.sP.fund_tc = value_float;
			}
			//Location of sites//////////////////////////////////////
			else if (lines.find("site_location") != std::string::npos) {
				sitesPos.insert(sitesPos.end(),value_vector.begin(),value_vector.end());
				currentSite++;

				//Adds variable mu for each sites
				int startB = lines.find('{', endP);
				int endB = lines.find('}', startB);
				if (startB > -1 || endB > -1) {
					double tempo = std::stof(lines.substr(startB + 1, endB - (startB + 1)));
					mu_array.push_back(tempo);
				}
			}
			//Allowed Jumps//////////////////////////////////////////////
			else if (lines.find("allow_jump") != std::string::npos) {
				allowedJump.insert(allowedJump.end(),value_vector.begin(),value_vector.end());
				numberOfAllowedJump++;
        
                //Adds the energy according to the jump
				int startB = lines.find('{', endP);
				int endB = lines.find('}', startB);
				double tempo = std::stof(lines.substr(startB + 1, endB - (startB + 1)));
                jumpEnergy.push_back(tempo);
			}
			//Get initial s states
			else if (lines.find("initial_states") != std::string::npos){
				std::vector<int> comaIndex;

				//Looks for all the commas and brakets
				comaIndex.push_back(startInfo);
				int i = startInfo;
				while(i>0){
					i = lines.find(',',i+1);
					if(i<0){
						break;
					}
					comaIndex.push_back(i);
					
				}
				comaIndex.push_back(endInfo);

				//What if there are no commas
				if(comaIndex.size()>2 || ((endInfo - startInfo)>1)){
					useSpecifiedR =true;
				}
				else{
					continue;
				}
				
				//Extract states
				std::vector<sType> tempInitialStateNum;
				for(unsigned int i = 0; i<comaIndex.size()-1;i++){
					std::string a = lines.substr(comaIndex.at(i) + 1, comaIndex.at(i+1) - (comaIndex.at(i) + 1));
					tempInitialStateNum.push_back(std::strtoll(a.c_str(),NULL,0));

				}	
				initState = tempInitialStateNum;
			}
			//Get initial k states
			else if (lines.find("added_site") != std::string::npos) {
				sendInfo.gP.g_AddedSite = value_int;
			}
			else if (lines.find("added_spin") != std::string::npos) {
				sendInfo.gP.g_AddedSpin = value_int;
			}
			else if (lines.find("eta") != std::string::npos) {
				sendInfo.gP.g_EtaValue = value_float;
			}
			else if (lines.find("add_all") != std::string::npos) {
				sendInfo.gP.g_AddAll = (value_str=="YES");
			}
			else if (lines.find("average") != std::string::npos) {
				sendInfo.gP.g_average = (value_str=="YES");
			}
			else if (lines.find("compute_green") != std::string::npos) {
				sendInfo.gP.g_compute = (value_str=="YES");
			}
		}
		myFile.close();
	}
	else { std::cout << "Unable to open the parameter file named : "<< fileRead<<std::endl; exit(1);}
	
	//Verification of entered parameters
	//n_sites = 0
	if(sendInfo.hubP.n_sites == 0) {std::cout<<"ERROR : The number of sites can't be zero"<<std::endl; exit(1);}
	// N > 2*sites
	if (N > sendInfo.hubP.n_sites * 2) { std::cout<<"ERROR : Too many electrons\nYou have put more electrons than twice the number of sites"<<std::endl; exit(1);}
	// sites != site position entered
	if(sendInfo.hubP.n_sites!=currentSite){
		std::cout<<"ERROR : The number of sites entered doesn't match the number of sites location"<<std::endl;
		exit(1);
	}
	//Percentage not between 0 and 1
	if (percentageOfStates < 0 || percentageOfStates > 1) {
		std::cout<<"ERROR : Percentage of states R entered not between 0 and 1"<<std::endl;
		exit(1);
	}
	
	float up_f = (N+Sz)/2;
	float down_f = (N-Sz)/2;
	if ((up_f != (int)up_f) || (down_f != (int)down_f)){
		std::cout<<"N and Sz are incompatible!"<<std::endl;	
		exit(1);
	}
	elec.up = up_f;
	elec.down = down_f;
	
	//Creating tMatrix
	double* tMatrix = new double[sendInfo.hubP.n_sites * sendInfo.hubP.n_sites]();
	for (unsigned char i = 0; i < sendInfo.hubP.n_sites; i++){	
		for (unsigned char j = i; j < sendInfo.hubP.n_sites; j++){
			if(i==j){
				tMatrix[i*sendInfo.hubP.n_sites+j]=0;
				continue;
			}
			int jumpX = sitesPos.at(i*3+0) - sitesPos.at(j*3+0); 	
			int jumpY = sitesPos.at(i*3+1) - sitesPos.at(j*3+1); 
			int jumpZ = sitesPos.at(i*3+2) - sitesPos.at(j*3+2);

			for(int k = 0; k < numberOfAllowedJump;k++){
				bool goodJump =	(jumpX == allowedJump.at(k*3+0) && jumpY == allowedJump.at(k*3+1) && jumpZ == allowedJump.at(k*3+2))\
								|| (jumpX == -allowedJump.at(k*3+0) && jumpY == -allowedJump.at(k*3+1) && jumpZ == -allowedJump.at(k*3+2));
				if(goodJump){
					tMatrix[i*sendInfo.hubP.n_sites+j] = tMatrix[j*sendInfo.hubP.n_sites+i] = jumpEnergy.at(k);
				}
			}
		}
	}
	
	sendInfo.hubP.tMatrix = std::vector<double>(tMatrix,tMatrix+(sendInfo.hubP.n_sites*sendInfo.hubP.n_sites));
	delete[] tMatrix;

	//Store mu_array and nu_matrix;
	
	init = createAntiFerro(sendInfo.hubP.n_sites,up_f,down_f);

	std::vector<sType> temp;
	temp.push_back(init);
	
	//Using entered or calculated states
	//R
	if(percentageOfStates > 1) percentageOfStates = 1;
	if (useSpecifiedR) sendInfo.sP.initState = initState;
	else sendInfo.sP.initState = temp;
	
	
	sendInfo.hubP.nElec = N;
	sendInfo.hubP.totSpin = Sz;

	if(sendInfo.sP.samplingSize <=0) sendInfo.sP.samplingSize = 1;

	return sendInfo;
	
}
