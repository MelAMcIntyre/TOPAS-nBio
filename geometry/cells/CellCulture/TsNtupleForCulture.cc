// Scorer for TsNtupleForCulture
//
// ********************************************************************
// *                                                                  *
// * This file is part of the TOPAS-nBio extensions to the            *
// *   TOPAS Simulation Toolkit.                                      *
// * The TOPAS-nBio extensions are freely available under the license *
// *   agreement set forth at: https://topas-nbio.readthedocs.io/     *
// *                                                                  *
// ********************************************************************
//
// Tuple scorer for scoring events deposited in the Culture organelles

#include "TsNtupleForCulture.hh"

#include "G4SystemOfUnits.hh"

#include <iomanip>
using std::setprecision;
#include <iostream>
#include <fstream>
#include <string>

TsNtupleForCulture::TsNtupleForCulture(TsParameterManager *pM, TsMaterialManager *mM, TsGeometryManager *gM, TsScoringManager *scM, TsExtensionManager *eM,
                                       G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
    : TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{

    // SetScorer();
    fNtuple->RegisterColumnF(&fPosX, "Position X", "cm");
    fNtuple->RegisterColumnF(&fPosY, "Position Y", "cm");
    fNtuple->RegisterColumnF(&fPosZ, "Position Z", "cm");
    fNtuple->RegisterColumnF(&fEnergy, "Energy", "MeV");
    fNtuple->RegisterColumnF(&fEnergyDep, "Energy Deposited", "MeV");
    fNtuple->RegisterColumnI(&fParticleType, "Particle Type (in PDG Format)");
    fNtuple->RegisterColumnI(&fTrackID, "Track ID");
    fNtuple->RegisterColumnI(&fRunID, "Run ID");
    fNtuple->RegisterColumnI(&fEventID, "Event ID");
    fNtuple->RegisterColumnD(&fGTime, "Global Time", "ps");
    fNtuple->RegisterColumnS(&fVolName, "Volume Name");
    fNtuple->RegisterColumnS(&fProcessName, "Process Name");

    std::ofstream rmfile;
    rmfile.open("output.csv");
    std::remove("output.csv");
    rmfile.close();

    std::ofstream rmeDepfile;
    rmeDepfile.open("energyDep.csv");
    std::remove("energyDep.csv");
    rmeDepfile.close();
}

TsNtupleForCulture::~TsNtupleForCulture() { ; }

G4bool TsNtupleForCulture::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
    ResolveSolid(aStep);

    G4double flagEnergyDep = aStep->GetTotalEnergyDeposit();

    G4StepPoint *theStepPoint = aStep->GetPreStepPoint();

    // Find volume name
    G4Track *aTrack = aStep->GetTrack();
    // G4String volumeName = aTrack->GetVolume()->GetName();
    G4VPhysicalVolume *volume = theStepPoint->GetTouchableHandle()->GetVolume();
    G4String volumeName = volume->GetName();

    // Get position
    G4ThreeVector pos = theStepPoint->GetPosition();

    // Score events that deposit energy in...
    G4String cylString = "Cyl";
    G4String cellString = "Cell";
    G4String nucString = "Nucleus";

    // Open output files
    std::ofstream outfile;
    outfile.open("output.csv", std::ofstream::app);
    std::ofstream eDepOutfile;
    eDepOutfile.open("energyDep.csv", std::ofstream::app);

    if ((flagEnergyDep >= 1.079e-5) && volumeName.contains(cylString))
    {
        // Get position
        fPosX = pos.x();
        fPosY = pos.y();
        fPosZ = pos.z();

        // Get copy numbers
        G4TouchableHistory *theTouchable = (G4TouchableHistory *)(theStepPoint->GetTouchable());
        fcopyNo = theTouchable->GetVolume()->GetCopyNo();        // Copy of Cylinder number
        fmotherCopyNo = theTouchable->GetVolume(2)->GetCopyNo(); // Copy of cell number

        // Get particle Energy
        fEnergy = theStepPoint->GetKineticEnergy();

        // Get Edep
        fEnergyDep = flagEnergyDep;

        // Get particle type
        fParticleType = aStep->GetTrack()->GetDefinition()->GetPDGEncoding();

        // Get IDs
        fTrackID = aStep->GetTrack()->GetTrackID();
        fRunID = GetRunID();
        fEventID = GetEventID();

        // Time
        fGTime = aStep->GetPreStepPoint()->GetGlobalTime();

        // Get volume Name
        fVolName = volumeName;

        // Get process name
        fProcessName = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
        G4double flagProcess = 0;

        if (fProcessName == "e-_G4DNAElastic")
            flagProcess = 11;
        else if (fProcessName == "e-_G4DNAExcitation")
            flagProcess = 12;
        else if (fProcessName == "e-_G4DNAIonisation")
            flagProcess = 13;
        else if (fProcessName == "e-_G4DNAAttachment")
            flagProcess = 14;
        else if (fProcessName == "e-_G4DNAVibExcitation")
            flagProcess = 15;
        else if (fProcessName == "eCapture")
            flagProcess = 16;
        else if (fProcessName == "proton_G4DNAExcitation")
            flagProcess = 21;
        else if (fProcessName == "proton_G4DNAIonisation")
            flagProcess = 22;
        else if (fProcessName == "proton_G4DNAChargeDecrease")
            flagProcess = 23;
        else if (fProcessName == "hydrogen_G4DNAExcitation")
            flagProcess = 31;
        else if (fProcessName == "hydrogen_G4DNAIonisation")
            flagProcess = 32;
        else if (fProcessName == "hydrogen_G4DNAChargeIncrease")
            flagProcess = 33;
        else if (fProcessName == "alpha_G4DNAExcitation")
            flagProcess = 41;
        else if (fProcessName == "alpha_G4DNAIonisation")
            flagProcess = 42;
        else if (fProcessName == "alpha_G4DNAChargeDecrease")
            flagProcess = 43;
        else if (fProcessName == "alpha+_G4DNAExcitation")
            flagProcess = 51;
        else if (fProcessName == "alpha+_G4DNAIonisation")
            flagProcess = 52;
        else if (fProcessName == "alpha+_G4DNAChargeDecrease")
            flagProcess = 53;
        else if (fProcessName == "alpha+_G4DNAChargeIncrease")
            flagProcess = 54;
        else if (fProcessName == "helium_G4DNAExcitation")
            flagProcess = 61;
        else if (fProcessName == "helium_G4DNAIonisation")
            flagProcess = 62;
        else if (fProcessName == "helium_G4DNAChargeIncrease")
            flagProcess = 63;

        outfile << fEventID << "," << fTrackID << "," << fmotherCopyNo << "," << fcopyNo << "," << fPosX << "," << fPosY << "," << fPosZ << "," << std::fixed << std::setprecision(16) << fGTime << "," << std::fixed << std::setprecision(1) << fParticleType << "," << flagProcess << "," << std::fixed << std::setprecision(6) << fEnergy << "," << fEnergyDep << "\n";
        fNtuple->Fill();
        return true;
    }

    if ((flagEnergyDep > 0) && volumeName.contains(cylString))
    {
        G4TouchableHistory *theTouchable = (G4TouchableHistory *)(theStepPoint->GetTouchable());
        fmotherCopyNo = theTouchable->GetVolume(2)->GetCopyNo(); // Copy of cell number
        eDepOutfile << fmotherCopyNo << "," << fEnergyDep << "\n";
    }
    else if ((flagEnergyDep > 0) && volumeName.contains(nucString))
    {
        // Get Edep
        fEnergyDep = flagEnergyDep;

        G4TouchableHistory *theTouchable = (G4TouchableHistory *)(theStepPoint->GetTouchable());
        fmotherCopyNo = theTouchable->GetVolume(1)->GetCopyNo(); // Copy of cell number

        eDepOutfile << fmotherCopyNo << "," << fEnergyDep << "\n";
    }
    else if ((flagEnergyDep > 0) && volumeName.contains(cellString))
    {
        fEnergyDep = flagEnergyDep;

        G4TouchableHistory *theTouchable = (G4TouchableHistory *)(theStepPoint->GetTouchable());
        fmotherCopyNo = theTouchable->GetVolume(0)->GetCopyNo(); // Copy of cell number

        eDepOutfile << fmotherCopyNo << "," << fEnergyDep << "\n";
    }
    outfile.close();
    eDepOutfile.close();

    return false;
}
