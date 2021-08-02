/*
 * man.cc
 *
 *  Created on: 18 Mar 2014
 *      Author: melda
 */

#include <QGraphicsWidget>
#include <QtGui>
#include <QtWidgets>
#include "MainWindow.h"
#include "GLWidget.h"

#include "../TissueFolding/SourceCode/Simulation.h"
#include <vector>

class MainWindow;
class SurfaceBase;
class GLWidget;

Simulation* Sim01;
int main(int argc, char **argv)
{
	bool Success = false;
	Sim01 = new Simulation();
	Sim01->displayIsOn = true;
	if (argc<2){
		Sim01->DisplaySave = false;
		cerr<<"Using default settings"<<endl;
		Success = true;
	}
	else{
		Success = Sim01->readExecutableInputs(argc, argv);
	}
	if (Success == 0 ){
		cout<<"Error in input to executable"<<endl;
		return true;
	}

	QApplication app(argc, argv);

	if (Sim01->DisplaySave){
		cout<<"Initiating simulation display"<<endl;
		Success = Sim01->initiateSavedSystem();
	}
	else{
		Success = Sim01->initiateSystem();
		if (Success == 0 ){
			cout<<"System is not initiated successfully, terminating"<<endl;
			return true;
		}
		else{
			cout<<"system initiated"<<endl;
		}
		int n = Sim01->Elements.size();
		for (int i=0; i<n; ++i){
			Sim01->Elements[i]->updatePositions(Sim01->Nodes);
		}
	}

	if (Success == 0 ){
		cout<<"System is not initiated successfully, terminating"<<endl;
		return true;
	}

	MainWindow mw(Sim01);
	mw.show();
	mw.MainGLWidget->show();
	mw.raise();
    //mw.setGeometry(1600, 600, 1600, 900); //ubuntu two screen
	//mw.setGeometry(1600, 600, 2000, 900); //ubuntu two screen
	mw.setGeometry(1600, 900, 1500, 900); //ubuntu two screen/no control panel

	//mw.setGeometry(50, 50, 900, 550); //mac
    mw.MainScene->update();

	return app.exec();
}




