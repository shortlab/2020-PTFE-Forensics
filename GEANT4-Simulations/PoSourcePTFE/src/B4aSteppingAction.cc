//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4aSteppingAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
// 
/// \file B4aSteppingAction.cc
/// \brief Implementation of the B4aSteppingAction class 

#include "B4aSteppingAction.hh"
#include "B4aEventAction.hh"
#include "B4DetectorConstruction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VParticleChange.hh"
#include "B4RunAction.hh"
#include "G4Decay.hh"
#include "G4RadioActiveDecay.hh"
#include "G4Run.hh"


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

G4Decay* theDecayProcess = new G4Decay();

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::B4aSteppingAction(
                      B4DetectorConstruction* detectorConstruction,
                      B4aEventAction* eventAction, B4RunAction* runAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction),
    fRunAction(runAction)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4aSteppingAction::~B4aSteppingAction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4aSteppingAction::UserSteppingAction(const G4Step* theStep)
{


  //get event #
  G4int eID = 0;
  const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
  const G4Run* run = G4RunManager::GetRunManager()->GetCurrentRun();
  G4int rID = run->GetRunID();
  if(evt) eID = evt->GetEventID();
  G4Track *theTrack = theStep->GetTrack(); 



  G4StepPoint* preStepPoint = theStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = theStep->GetPostStepPoint();

// Get incident position and energy
if (theTrack->GetParticleDefinition()->GetParticleName() == "alpha" & postStepPoint->GetProcessDefinedStep() != 0) {
          if (preStepPoint->GetPhysicalVolume()->GetName() == "World" & postStepPoint->GetPhysicalVolume()->GetName() == "0") {
                G4double xposition_um = (postStepPoint->GetPosition().getX()/um);
                G4double yposition_um = (postStepPoint->GetPosition().getY()/um);  
                G4double energy_MeV = postStepPoint->GetKineticEnergy()/MeV;
                G4String fileName = "Position_Incident_PTFE_um.txt";
                std::ostringstream commandOS;
                commandOS << fileName;
                std::ofstream ofile;
                ofile.open (G4String(commandOS.str()), ios::out | ios::app);     // ascii file   
                // in cm and MeV 
                ofile << eID << " " << theTrack->GetParticleDefinition()->GetParticleName() << " " << energy_MeV << " " << xposition_um << " " << yposition_um << "\n";
                ofile.close(); 
               // theTrack->SetTrackStatus(fKillTrackAndSecondaries); 
              } } 


// Get energy deposited in depth bins
if (postStepPoint->GetProcessDefinedStep() != 0 & theTrack->GetVolume()->GetName() != "World" & theTrack->GetVolume()->GetName() != "SourceLayer" & theTrack->GetVolume()->GetName() != "AuLayer") {
                G4double xposition_um = (postStepPoint->GetPosition().getX()/um);
                G4double yposition_um = (postStepPoint->GetPosition().getY()/um);  
                G4double zposition_um = (postStepPoint->GetPosition().getZ()/um);  

                G4double eDep_MeV = theStep->GetTotalEnergyDeposit()/MeV;
                G4String fileName = "Edep_PTFE_MeV_um.txt";
                std::ostringstream commandOS;
                commandOS << fileName;
                std::ofstream ofile;
                ofile.open (G4String(commandOS.str()), ios::out | ios::app);     // ascii file   
                // in cm and MeV 
                ofile << eID << " " << theTrack->GetParticleDefinition()->GetParticleName() << " " << eDep_MeV << " " << xposition_um << " " << yposition_um << " " << theTrack->GetVolume()->GetName() << " \n";
                ofile.close(); 
     } 






}























//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......