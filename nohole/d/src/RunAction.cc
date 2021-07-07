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
/// \file optical/OpNovice2/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Make this appear first!
#include "G4Timer.hh"

#include "RunAction.hh"
//#include "RunData.hh"
#include "HistoManager.hh"
#include "PrimaryGeneratorAction.hh"

#include "Run.hh"
#include "G4Run.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(PrimaryGeneratorAction* prim)
 : G4UserRunAction(),
   fTimer(nullptr),
   fRun(nullptr),
   fHistoManager(nullptr),
   fPrimary(prim)
{
  fTimer = new G4Timer;


  auto  analysisManager = G4AnalysisManager::Instance();
  // analysisManager->CreateH1("sipm1_photons","Counts", 500, 0., 500);
  // analysisManager->CreateH1("sipm1_energy","Counts", 10, 0., 10*eV);
  // analysisManager->CreateH1("sipm2_photons","Counts", 500, 0., 500);
   //analysisManager->CreateH1("sipm2_energy","Counts", 10, 0., 10*eV);
   // analysisManager->CreateH1("Edep","count", 300, 0., 800*keV); 

  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);

  analysisManager->CreateNtuple("dosimetry", "");

  for (G4int i_sipm = 0; i_sipm < 2; i_sipm++) {
        buffer = "Count" + to_string(i_sipm + 1);
    analysisManager->CreateNtupleDColumn(buffer);
    buffer = "Ekin_eV" + to_string(i_sipm + 1);
    analysisManager->CreateNtupleDColumn(buffer);
    }

    analysisManager->CreateNtupleDColumn("EnterancexAxes");
    analysisManager->CreateNtupleDColumn("EnteranceyAxes");
    analysisManager->CreateNtupleDColumn("EnterancezzAxes");

    analysisManager->CreateNtupleDColumn("ExitxAxes");
    analysisManager->CreateNtupleDColumn("ExityAxes");
    analysisManager->CreateNtupleDColumn("ExitzAxes");

    analysisManager->CreateNtupleDColumn("ScixAxes");
    analysisManager->CreateNtupleDColumn("SciyAxes");
    analysisManager->CreateNtupleDColumn("ScizAxes");

    analysisManager->CreateNtupleDColumn("EnteranceTheta");
    analysisManager->CreateNtupleDColumn("EnterancePhi");
        analysisManager->CreateNtupleDColumn("ExitTheta");
        analysisManager->CreateNtupleDColumn("ExitPhi");
        //
analysisManager->CreateNtupleDColumn("TotalEnergyDeposit_keV");
analysisManager->CreateNtupleDColumn("Dose_Gy");
analysisManager->CreateNtupleDColumn("Dose_Equivalent_Sv");
//

analysisManager->CreateNtupleDColumn("Incident_electron_spectrum_keV");
analysisManager->CreateNtupleDColumn("Incident_gamma_spectrum_keV");
analysisManager->CreateNtupleDColumn("Incident_alpha_spectrum_keV");
analysisManager->CreateNtupleDColumn("Incident_proton_spectrum_keV");
analysisManager->CreateNtupleDColumn("Incident_neutron_spectrum_keV");
analysisManager->CreateNtupleDColumn("Incident_tot_spectrum_keV");
//
analysisManager->CreateNtupleDColumn("Ent_electron_spectrum_keV");
analysisManager->CreateNtupleDColumn("Ent_gamma_spectrum_keV");
analysisManager->CreateNtupleDColumn("Ent_alpha_spectrum_keV");
analysisManager->CreateNtupleDColumn("Ent_proton_spectrum_keV");
analysisManager->CreateNtupleDColumn("Ent_neutron_spectrum_keV");
analysisManager->CreateNtupleDColumn("Ent_tot_spectrum_keV");
//
analysisManager->CreateNtupleDColumn("Exit_electron_spectrum_keV");
analysisManager->CreateNtupleDColumn("Exit_gamma_spectrum_keV");
analysisManager->CreateNtupleDColumn("Exit_alpha_spectrum_keV");
analysisManager->CreateNtupleDColumn("Exit_proton_spectrum_keV");
analysisManager->CreateNtupleDColumn("Exit_neutron_spectrum_keV");
analysisManager->CreateNtupleDColumn("Exit_tot_spectrum_keV");
//
analysisManager->CreateNtupleDColumn("Incident_electron_count");
analysisManager->CreateNtupleDColumn("Incident_gamma_count");
analysisManager->CreateNtupleDColumn("Incident_neutron_count");
analysisManager->CreateNtupleDColumn("Incident_proton_count");
analysisManager->CreateNtupleDColumn("Incident_alpha_count");
//
analysisManager->CreateNtupleDColumn("Ent_electron_count");
analysisManager->CreateNtupleDColumn("Ent_gamma_count");
analysisManager->CreateNtupleDColumn("Ent_neutron_count");
analysisManager->CreateNtupleDColumn("Ent_proton_count");
analysisManager->CreateNtupleDColumn("Ent_alpha_count");
//
analysisManager->CreateNtupleDColumn("Exit_electron_count");
analysisManager->CreateNtupleDColumn("Exit_gamma_count");
analysisManager->CreateNtupleDColumn("Exit_neutron_count");
analysisManager->CreateNtupleDColumn("Exit_proton_count");
analysisManager->CreateNtupleDColumn("Exit_alpha_count");
analysisManager->CreateNtupleDColumn("Photons_count");
//
analysisManager->CreateNtupleDColumn("coincidence_photon");
analysisManager->CreateNtupleDColumn("coincidence_energy_keV");

  analysisManager->FinishNtuple();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fTimer;
  delete G4AnalysisManager::Instance();
}

G4Run* RunAction::GenerateRun()
{
  fRun = new Run();
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  if (fPrimary) {
    G4ParticleDefinition* particle =
      fPrimary->GetParticleGun()->GetParticleDefinition();
    G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
    fRun->SetPrimary(particle, energy);
  }

    auto analysisManager = G4AnalysisManager::Instance();
    G4String fileName = "dosimetry";
    analysisManager->OpenFile(fileName);

    fTimer->Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  fTimer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent()
         << " " << *fTimer << G4endl;

  if (isMaster) fRun->EndOfRun();

  // save histograms
    auto analysisManager = G4AnalysisManager::Instance();

    analysisManager->Write();
    analysisManager->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
