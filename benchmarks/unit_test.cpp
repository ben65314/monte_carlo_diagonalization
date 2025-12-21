#include "test_basicFunctions.h"
#include "test_electronManipulationFunctions.h"
//#include "test_advancedFunctions.h"
//#include "test_lanczosAlgorithm.h"
//#include "test_statesArr.h"

int main(){
	bool details = true;
	allTests_basicFunctions(details);
	allTests_electronManipulationFunctions(details);
	//allTests_advancedFunctions(details);
	//allTests_lanczosAlgorithm(details);
	//allTest_StatesArr_for_childs(details);
    return 0;
}
