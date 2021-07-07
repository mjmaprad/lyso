/// \file RunData.cc
/// \brief Implementation of the RunData class

#include "RunData.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EventManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunData::RunData()
        : G4Run()
{

fposx_sci=-9;
     fposy_sci=-9;
     fposz_sci=-9;
     fposx_ent=-9;
     fposy_ent=-9;
     fposz_ent=-9;
     fposx_exi=-9;
     fposy_exi=-9;
     fposz_exi=-9;
     phi_ent=-9;
     theta_ent=-9;
     phi_exi=-9;
     theta_exi=-9;
  

     eDep=0;
     dose=0;
    eq_dose=0;
     fcount_photon=0;
     ////
     elecspec=0;
     gammaspec=0;
     protonspec=0;
 neutronspec=0;
     alphaspec=0;
     totspec=0;    //
     ////////
          entelecspec=0;
     entgammaspec=0;
     entprotonspec=0;
 entneutronspec=0;
     entalphaspec=0;
     enttotspec=0;    //
         entelecnum=0;
     entgammanum=0;
     entneutronnum=0;
     entalphanum=0;
          entprotonnum=0;
     /////////////////
          ////////
          extelecspec=0;
     extgammaspec=0;
     extprotonspec=0;
 extneutronspec=0;
     extalphaspec=0;
     exttotspec=0;    //
         extelecnum=0;
     extgammanum=0;
     extneutronnum=0;
     extalphanum=0;
          extprotonnum=0;
     /////////////////
      elecnum=0;
     gammanum=0;
     neutronnum=0;
     alphanum=0;
          protonnum=0;
for (G4int i_sipm = 0; i_sipm < 2; i_sipm++) {
  fen[i_sipm] = 0;
  fcount[i_sipm] = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunData::~RunData()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunData::FillPerEvent()
{
//   G4cout
// <<G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()
// <<G4endl;
    G4int counter = 0;

 auto analysisManager = G4AnalysisManager::Instance();

   // analysisManager->FillH1(0, fcount[0]);
    //analysisManager->FillH1(1, fen[0]);   
   // analysisManager->FillH1(1, fcount[1]);
    // analysisManager->FillH1(3, fen[1]); 
    // analysisManager->FillH1(4, eDep);
    // accumulate statistic
    // in the order od the histograms, ntuple columns declarations
    for (G4int i_sipm = 0; i_sipm < 2; i_sipm++) {
   analysisManager->FillNtupleDColumn(counter++,fcount[i_sipm]);
    analysisManager->FillNtupleDColumn(counter++,fen[i_sipm]);

    }
    analysisManager->FillNtupleDColumn(counter++,fposx_ent);
    analysisManager->FillNtupleDColumn(counter++,fposy_ent);
    analysisManager->FillNtupleDColumn(counter++,fposz_ent);
    ////
    analysisManager->FillNtupleDColumn(counter++,fposx_exi);
    analysisManager->FillNtupleDColumn(counter++,fposy_exi);
    analysisManager->FillNtupleDColumn(counter++,fposz_exi);
    //
    analysisManager->FillNtupleDColumn(counter++,fposx_sci);
    analysisManager->FillNtupleDColumn(counter++,fposy_sci);
    analysisManager->FillNtupleDColumn(counter++,fposz_sci);
    ///
        analysisManager->FillNtupleDColumn(counter++,theta_ent);
    analysisManager->FillNtupleDColumn(counter++,phi_ent);
        analysisManager->FillNtupleDColumn(counter++,theta_exi);
        analysisManager->FillNtupleDColumn(counter++,phi_exi);
        ////
    analysisManager->FillNtupleDColumn(counter++,eDep);
        analysisManager->FillNtupleDColumn(counter++,dose);
              analysisManager->FillNtupleDColumn(counter++,eq_dose);
//
    analysisManager->FillNtupleDColumn(counter++,elecspec);
    analysisManager->FillNtupleDColumn(counter++,gammaspec);
            analysisManager->FillNtupleDColumn(counter++,alphaspec);
                    analysisManager->FillNtupleDColumn(counter++,protonspec);
                analysisManager->FillNtupleDColumn(counter++,neutronspec);
                analysisManager->FillNtupleDColumn(counter++,totspec);
//
 analysisManager->FillNtupleDColumn(counter++,entelecspec);
    analysisManager->FillNtupleDColumn(counter++,entgammaspec);
            analysisManager->FillNtupleDColumn(counter++,entalphaspec);
                    analysisManager->FillNtupleDColumn(counter++,entprotonspec);
                analysisManager->FillNtupleDColumn(counter++,entneutronspec);
                analysisManager->FillNtupleDColumn(counter++,enttotspec);
//                                
//
 analysisManager->FillNtupleDColumn(counter++,extelecspec);
    analysisManager->FillNtupleDColumn(counter++,extgammaspec);
            analysisManager->FillNtupleDColumn(counter++,extalphaspec);
                    analysisManager->FillNtupleDColumn(counter++,extprotonspec);
                analysisManager->FillNtupleDColumn(counter++,extneutronspec);
                analysisManager->FillNtupleDColumn(counter++,exttotspec);
 //

                //
                  analysisManager->FillNtupleDColumn(counter++,elecnum);
                 analysisManager->FillNtupleDColumn(counter++,gammanum);
               analysisManager->FillNtupleDColumn(counter++,neutronnum);
                analysisManager->FillNtupleDColumn(counter++,protonnum);
                analysisManager->FillNtupleDColumn(counter++,alphanum);
 //

                //
                  analysisManager->FillNtupleDColumn(counter++,entelecnum);
                 analysisManager->FillNtupleDColumn(counter++,entgammanum);
               analysisManager->FillNtupleDColumn(counter++,entneutronnum);
                analysisManager->FillNtupleDColumn(counter++,entprotonnum);
                analysisManager->FillNtupleDColumn(counter++,entalphanum);
 //

                //
                  analysisManager->FillNtupleDColumn(counter++,extelecnum);
                 analysisManager->FillNtupleDColumn(counter++,extgammanum);
               analysisManager->FillNtupleDColumn(counter++,extneutronnum);
                analysisManager->FillNtupleDColumn(counter++,extprotonnum);
                analysisManager->FillNtupleDColumn(counter++,extalphanum);
                analysisManager->FillNtupleDColumn(counter++,fcount_photon);


            if (fcount[0]!=0&&fcount[1]!=0)  {
 analysisManager->FillNtupleDColumn(counter++,fcount[0]+fcount[1]);
  analysisManager->FillNtupleDColumn(counter++,(fen[0]+fen[1])/1000);
            }
             else  analysisManager->FillNtupleDColumn(counter++,0);

    analysisManager->AddNtupleRow();
         
  Reset();




// //G4cout<<num_Cphotons_ck<<"  "<<num_Cphotons<<G4endl;

//   if ((G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()+1)==totNumEvent())
 
//   {
//      G4cout<<"- Total dose(Gy)    Dose equivalent rate (Sv)  Dose equivalent rate (Sv/h)-----------"<<G4endl;
//       G4cout<<accdose<<" "<<acc_eq_dose<<" "<<acc_eq_dose*3600<<G4endl;
//        G4cout<<"---------------------------------total photons created-------------------------------"<<G4endl;
//      G4cout<<accphoton<<G4endl;
//        G4cout<<"---------------------------------entered particles-------------------------------"<<G4endl;
//     //  G4cout<<tot_fcoince<<G4endl;
//         G4cout<<"- No.of incident e gamma alpha neutron proton------------------------------"<<G4endl;
//       G4cout<<accelecnum<<" "<<accgammanum<<" "<<accalphanum<<" "<<accneutronnum<<" "<<accprotonnum<<G4endl;
//         G4cout<<"---------------------------------total energy in sipm---------------------"<<G4endl;

//         G4cout<<"In sipm1: "<<tot_sipme1/1e3<<" keV"<<G4endl;
//         G4cout<<"In sipm2: "<<tot_sipme2/1e3<<" keV"<<G4endl;
//        G4cout<<"---------------------------------cpm in sipm-------------------------------"<<G4endl;
//         G4cout<<"In sipm1: "<<tot_sipmc1<<G4endl;
//         G4cout<<"In sipm2: "<<tot_sipmc2<<G4endl;
      
//   }
//tot_fcoince ++;
//  if (G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()==0){
//       tot_dose=0;
//     tot_photon=0;
//     tot_eq_dose=0;
//          accelecnum=0;
//      accgammanum=0;
//      accneutronnum=0;
//      accalphanum=0;
//           accprotonnum=0;
//                tot_sipmcount1=0;
//      tot_sipmenergy1=0;
//        tot_sipmcount2=0;
//     tot_sipmenergy2=0;
//  }
    


// G4cout<<tot_dose<<" "<<tot_photon<<" "<<tot_eq_dose<<""<<G4endl;
//     //tot_dose+=dose;
//   //  tot_eq_dose+=eq_dose;
// //tot_photon+=fcount_photon;
//   //  accelecnum+=elecnum;
//    //  accprotonnum+=protonnum;
//     //accgammanum+=gammanum;
//     //accneutronnum+=neutronnum;
//     //accalphanum+=alphanum;
//  analysisManager->AddNtupleRow();

//  if ((G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID()+1)==
//         G4Run::GetNumberOfEventToBeProcessed())
//  {
//      G4cout<<"- Total dose(Gy)    Dose equivalent rate (Sv)  Dose equivalent rate (Sv/h)-----------"<<G4endl;
//       G4cout<<tot_dose<<" "<<tot_eq_dose<<" "<<tot_eq_dose*3600<<G4endl;
//        G4cout<<"---------------------------------total photons created-------------------------------"<<G4endl;
//      G4cout<<tot_photon<<G4endl;
//        G4cout<<"---------------------------------coincidence coefficient in the events-------------------------------"<<G4endl;
//      G4cout<<tot_fcoince<<G4endl;
//         G4cout<<"- No.of incident e gamma alpha neutron proton------------------------------"<<G4endl;
//       G4cout<<accelecnum<<" "<<accgammanum<<" "<<accalphanum<<" "<<accneutronnum<<" "<<accprotonnum<<G4endl;
//         G4cout<<"---------------------------------total energy in sipm---------------------"<<G4endl;

//         G4cout<<"In sipm1: "<<tot_sipmenergy1/1e3<<" keV"<<G4endl;
//         G4cout<<"In sipm2: "<<tot_sipmenergy2/1e3<<" keV"<<G4endl;
//        G4cout<<"---------------------------------cpm in sipm-------------------------------"<<G4endl;
//         G4cout<<"In sipm1: "<<tot_sipmcount1<<" keV"<<G4endl;
//         G4cout<<"In sipm2: "<<tot_sipmcount2<<" keV"<<G4endl;
      
      //  G4cout<<"- Total dose(Gy)   Dose equivalent rate (Sv/h)------------------------------"<<G4endl;
      // G4cout<<tot_dose<<" "<<tot_eq_dose*3600<<G4endl;
      //   G4cout<<"- No.of incident e gamma alpha neutron proton------------------------------"<<G4endl;
//      G4cout<<elecnum<<" "<<gammanum<<" "<<alphanum<<" "<<neutronnum<<" "<<protonnum<<G4endl;
//        G4cout<<"---------------------------------total energy in sipm---------------------"<<G4endl;
//            for (G4int i_sipm = 0; i_sipm <(sizeof(tot_sipm_energy)/sizeof(*tot_sipm_energy)); i_sipm++) {
//         G4cout<<"In sipm"<<i_sipm+1<<": "<<tot_sipm_energy[i_sipm]/1e3<<" keV"<<G4endl;}
     
//        G4cout<<"---------------------------------cpm in sipm-------------------------------"<<G4endl;
//            for (G4int i_sipm = 0; i_sipm <(sizeof(tot_sipm_count)/sizeof(*tot_sipm_count)); i_sipm++) {
//         G4cout<<"In sipm"<<i_sipm+1<<": "<<tot_sipm_count[i_sipm]*60<<" cpm"<<G4endl;}
//                G4cout<<"---------------------------------cpm in sipm-------------------------------"<<G4endl;
    //    G4cout<<"---------------------------------total photons created-------------------------------"<<G4endl;
    //  G4cout<<tot_photon<<G4endl;
    //    G4cout<<"---------------------------------coincidence coefficient in the events-------------------------------"<<G4endl;
    //  G4cout<<tot_fcoince<<G4endl;
 //}




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunData::Reset()
{
for (G4int i_sipm = 0; i_sipm < 2; i_sipm++) {
  fen[i_sipm] = 0;
  fcount[i_sipm] = 0;
  }
   fposx_sci=-9;
     fposy_sci=-9;
     fposz_sci=-9;
     fposx_ent=-9;
     fposy_ent=-9;
     fposz_ent=-9;
     fposx_exi=-9;
     fposy_exi=-9;
     fposz_exi=-9;
     phi_ent=-9;
     theta_ent=-9;
     phi_exi=-9;
     theta_exi=-9;
  
     eDep=0;
     dose=0;
    eq_dose=0;
     fcount_photon=0;
     ////
     elecspec=0;
     gammaspec=0;
     protonspec=0;
 neutronspec=0;
     alphaspec=0;
     totspec=0;    //
     ////////
          entelecspec=0;
     entgammaspec=0;
     entprotonspec=0;
 entneutronspec=0;
     entalphaspec=0;
     enttotspec=0;    //
         entelecnum=0;
     entgammanum=0;
     entneutronnum=0;
     entalphanum=0;
          entprotonnum=0;
     /////////////////
          ////////
          extelecspec=0;
     extgammaspec=0;
     extprotonspec=0;
 extneutronspec=0;
     extalphaspec=0;
     exttotspec=0;    //
         extelecnum=0;
     extgammanum=0;
     extneutronnum=0;
     extalphanum=0;
          extprotonnum=0;
     /////////////////
      elecnum=0;
     gammanum=0;
     neutronnum=0;
     alphanum=0;
          protonnum=0;
          for (G4int i_sipm = 0; i_sipm < 2; i_sipm++) {
  fen[i_sipm] = 0;
  fcount[i_sipm] = 0;
  }
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
