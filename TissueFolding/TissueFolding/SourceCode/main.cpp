using namespace std;

#include "Simulation.h"
#include <vector>

int main(int argc, char * argv[])
{
	Simulation* Sim01 = new Simulation();
	Sim01->displayIsOn = false;
	if (argc<2){
		Sim01->DisplaySave = false;
		cerr<<"Using default settings"<<endl;
	}
	else{
		bool Success = Sim01->readExecutableInputs(argc, argv);
		if (!Success){
			cerr<<"Error in input to executable"<<endl;
			return 1;
		}
	}
	if (Sim01->DisplaySave){
		cerr<<"This is the executable for running the model without display"<<endl;
		return true;
	}
	else{
		Sim01->initiateSystem();
		int n = Sim01->Elements.size();
		for (int i=0; i<n; ++i){
			Sim01->Elements[i]->updatePositions(Sim01->Nodes);
		}
		cout<<"Initiating simulation in the background"<<endl;
		while (Sim01->currSimTimeSec <= Sim01->SimLength){
			//cout<<"running step: "<<Sim01.timestep<<", this is time: "<<Sim01.timestep*Sim01.dt<<" sec"<<endl;
			bool Success = Sim01->runOneStep();
			if (Success == false ){
				break;
			}
		}
		Sim01->wrapUpAtTheEndOfSimulation();
		Sim01->writeRelaxedMeshFromCurrentState();
	}

	delete Sim01;
	cout<<"Finished Simulation"<<endl;
}

