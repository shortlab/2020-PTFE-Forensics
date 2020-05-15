#include "TrackingAction.hh"
#include "B4aEventAction.hh"
#include "B4aSteppingAction.hh"
#include "B4RunAction.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
//#include "AIDA/AIDA.h"



#include "B4aSteppingAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "B4RunAction.hh"
#include "G4Decay.hh"
#include "G4RadioActiveDecay.hh"


#include "G4ElectronIonPair.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

#include "Randomize.hh"
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdio.h> 
#include <math.h>


using namespace std;



TrackingAction::TrackingAction(B4aEventAction* EvAct, B4RunAction* run)
:evAction(EvAct), Run(run)
{ }


TrackingAction::~TrackingAction()
{ }

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
	const G4ParticleDefinition * particleDef = aTrack->GetParticleDefinition();
	G4String particleName = particleDef->GetParticleName();
	G4int atomicNumber = particleDef->GetAtomicNumber();


	/*if(atomicNumber > 1) {

        

          std::ostringstream commandOS;
          commandOS << "CaF2Window_Shielding_63MeV_10000000_pencil_Decay_WithTimes.txt";
          // 1 if true
          G4bool isRadioactive = !(aTrack->GetParticleDefinition()->GetPDGStable());
          // in seconds
          G4double lifetime = 0.69314718056*(aTrack->GetParticleDefinition()->GetPDGLifeTime()/second);
          G4double creationTime = aTrack->GetGlobalTime()/second;


          //if ((isRadioactive == true) && (lifetime > 1.0) && (lifetime < (60.*60.*24.*365.*100.))) {

          	//get event #
			  G4int eID = 0;
			  const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
			  const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
			  if(evt) eID = evt->GetEventID();
        
          std::ofstream ofile;
          ofile.open (G4String(commandOS.str()), ios::out | ios::app);     // ascii file    
          ofile << eID << " " << aTrack->GetTrackID() << " " << particleName << " " << aTrack->GetCreatorProcess()->GetProcessName() << " " << creationTime << " " << lifetime << " " << isRadioactive << endl;

          commandOS.str("");

          
        }
*/

}


void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{

}



