// Component for TsSphericalCell
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
// A simple spherical cell.
// User has the option of including organelles: nucleus and/or mitochondria.

#include "TsSphericalCell.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4LogicalVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "G4VSolid.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"

#include "globals.hh"
#include <math.h>

TsSphericalCell::TsSphericalCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM, TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name) :
TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}


TsSphericalCell::~TsSphericalCell()
{;}

void TsSphericalCell::ResolveParameters() {
    
    CellRadius = fPm->GetDoubleParameter(GetFullParmName("CellRadius"), "Length");    


}


G4VPhysicalVolume* TsSphericalCell::Construct()
{
	BeginConstruction();
    
    //***********************************************************************
    //              Envelope Geometry : spherical cell
    //***********************************************************************
    
    G4Orb* gCell = new G4Orb(fName, CellRadius);
    fEnvelopeLog = CreateLogicalVolume(gCell);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);
    
    
    //***********************************************************************
    // Optional : include a nucleus and/or mitochondria in the cell
    //***********************************************************************
    
    //***************************
    // Subcomponent: Nucleus
    //***************************
    
    G4double NuclRadius = 0.0*um;
    G4String name = GetFullParmName("Nucleus/NucleusRadius");
    if (fPm->ParameterExists(name)) {
        
        NuclRadius = fPm->GetDoubleParameter(name, "Length");
        G4String subComponentName1 = "Nucleus";
    
        G4RotationMatrix* rotNuc = new G4RotationMatrix();
        
        rotNuc->rotateX(0);
        rotNuc->rotateY(0);
        
        G4double transNucX = 0 * um;
        G4double transNucY = 0 * um;
        G4double transNucZ = 0 * um;
        
        G4String name1 = GetFullParmName("Nucleus/translateX");
        if (fPm -> ParameterExists(name1)){
            transNucX = fPm->GetDoubleParameter(name1, "Length");
            if (transNucX > NuclRadius) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                exit(1);
            }
        }
        
        name1 = GetFullParmName("Nucleus/translateY");
        if (fPm -> ParameterExists(name1)){
            transNucY = fPm->GetDoubleParameter(name1, "Length");
            if (transNucY > NuclRadius) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                exit(1);
            }
        }
        
        name1 = GetFullParmName("Nucleus/translateZ");
        if (fPm -> ParameterExists(name1)){
            transNucZ = fPm->GetDoubleParameter(name1, "Length");
            if (transNucZ > NuclRadius) {
                G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
                G4cerr << "Parameter " << name1 << " sets nucleus outside of cell." << G4endl;
                exit(1);
            }
            
        }
        
        
        G4ThreeVector* NucPos = new G4ThreeVector(transNucX,transNucY,transNucZ);
        
        G4Orb* gNucleus = new G4Orb("gNucleus", NuclRadius);
        G4LogicalVolume* lNucleus = CreateLogicalVolume(subComponentName1, gNucleus);
        G4VPhysicalVolume* pNucleus = CreatePhysicalVolume(subComponentName1, lNucleus, rotNuc, NucPos, fEnvelopePhys);
        
        G4bool OverlapCheck = pNucleus->CheckOverlaps();
        
        if (OverlapCheck == true){
            G4cerr << "Topas is exiting due to a serious error in geometry setup." << G4endl;
            G4cerr << "Nucleus overlaps with the cell." << G4endl;
            exit(1);
        }
        

        // fCreate = true;
        // if (fPm->ParameterExists(GetFullParmName("GenerateCylinders")))
        //         fCreate = fPm->GetBooleanParameter(GetFullParmName("GenerateCylinders"));
        
        // fOutputFile = "RandomCylinders_positions_nm_rotations_deg.txt";
        // if ( fPm->ParameterExists(GetFullParmName("OutputFile") ))
        //         fOutputFile = fPm->GetStringParameter(GetFullParmName("OutputFile"));

        // if ( !fCreate ) 
        //         fOutputFile = fPm->GetStringParameter(GetFullParmName("InputFile"));
        
        // fNumberOfCylinders = 10000;
        // if ( fPm->ParameterExists(GetFullParmName("NumberOfCylinders") ))
        //         fNumberOfCylinders = fPm->GetIntegerParameter(GetFullParmName("NumberOfCylinders"));
        
        // fTRMax = 0.5*2.3*nm;
        // fTHL   = 0.5*3.4*nm;
        // if (fPm->ParameterExists(GetFullParmName("Cylinders/RMax")))
        //         fTRMax = fPm->GetDoubleParameter(GetFullParmName("Cylinders/RMax"), "Length");
        // if ( fPm->ParameterExists(GetFullParmName("Cylinders/HL")))
        //         fTHL = fPm->GetDoubleParameter(GetFullParmName("Cylinders/HL"), "Length");
        
  
        // fRMax = fPm->GetDoubleParameter(GetFullParmName("RMax"), "Length");
        
        
        // if ( fCreate && fPm->GetIntegerParameter("Ts/NumberOfThreads") > 1 ) {
        //         G4cerr << "TOPAS is exiting due to an error. Conflict between component TsRandomCylindersInComponent "
        //         << "and Ts/NumberOfThreads. " << GetFullParmName("GenerateCylinders") << " is set to true. Then "
        //         << "Ts/NumberOfThredas must be set to 1." << G4endl;
        //         exit(1);
        // }

    //     G4double extensionX, extensionY, extensionZ;
        
    //     extensionX = 2.0 * fRMax;
    //     extensionY = extensionX;
    //     extensionZ = extensionY;

    //     G4String subComponentName2 = "Cylinders";
    //     G4Tubs* cylinder = new G4Tubs(subComponentName2, 0, fTRMax, fTHL, 0*deg, 360*deg);
    //     G4LogicalVolume* cylinderLog = CreateLogicalVolume(subComponentName2, cylinder);
        
    //     G4double x, y, z, u, v, w;

    //     if ( fCreate ) {
    //         std::ofstream outFile(fOutputFile);
    //         G4int n = 0;
    //         G4int volID = 1;
                
    //         while( n < fNumberOfCylinders ) {
    //             x = G4RandFlat::shoot(-0.5*extensionX,0.5*extensionX);
    //             y = G4RandFlat::shoot(-0.5*extensionY,0.5*extensionY);
    //             z = G4RandFlat::shoot(-0.5*extensionZ,0.5*extensionZ);
                
    //             EInside test_status = gNucleus->Inside(G4ThreeVector(x, y, z));

    //             if ( test_status == kInside ) {
    //                 G4RotationMatrix* rot = new G4RotationMatrix();
    //                 u = G4RandFlat::shoot(0., 180.0*deg);
    //                 v = G4RandFlat::shoot(0., 180.0*deg);
    //                 w = G4RandFlat::shoot(0., 180.0*deg);
    //                 rot->rotateX(u);
    //                 rot->rotateY(v);
    //                 rot->rotateZ(w);
    //                 G4ThreeVector* pos = new G4ThreeVector(x, y, z);
    //                 G4VPhysicalVolume* phys = CreatePhysicalVolume(subComponentName2, volID, true, cylinderLog,
    //                                                             rot, pos, pNucleus);
    //                 if ( phys->CheckOverlaps(1000, 0., false, 1) ) {
    //                     delete phys;
    //                 } else {
    //                     outFile << x/nm << " " << y/nm << " " << z/nm << " "
    //                             << u/deg << " " << v/deg << " " << w/deg << G4endl;
                        
    //                     volID++;
    //                     n++;
    //                 }
    //             }
    //         }
    //     } else {
    //         G4int volID = 1;
    //         std::ifstream inFile(fOutputFile);
    //         while(1) {
    //             inFile >> x >> y >> z >> u >> v >> w;
    //             if ( !inFile.good() ) break;
    //             G4RotationMatrix* rot = new G4RotationMatrix();
    //             rot->rotateX(u*deg);
    //             rot->rotateY(v*deg);
    //             rot->rotateZ(w*deg);
    //             CreatePhysicalVolume(subComponentName2, volID, true, cylinderLog, rot, new G4ThreeVector(x*nm,y*nm,z*nm), pNucleus);
    //             volID++;
    //         }
            
    //     }
    // }
    //*******************************
    // Subcomponent: Mitochondria
    //*******************************
    
    G4String name = GetFullParmName("Mitochondria/NumberOfMitochondria");
    if (fPm->ParameterExists(name)) {
        
        //number of mitochondria
        const G4int NbOfMito  = fPm->GetIntegerParameter( GetFullParmName("Mitochondria/NumberOfMitochondria") );
        
        //Semi-axis lengths of the ellpsoid/mitochondria (default values if none are specified)
        G4double EllA = 0.5*micrometer;
        G4double EllB = 0.3*micrometer;
        G4double EllC = 0.9*micrometer;
        
        name=GetFullParmName("Mitochondria/a");
        if (fPm->ParameterExists(name)){EllA = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/a"), "Length" );}
        
        name=GetFullParmName("Mitochondria/b");
        if (fPm->ParameterExists(name)){EllB = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/b"), "Length" );}
        
        name=GetFullParmName("Mitochondria/c");
        if (fPm->ParameterExists(name)){EllC = fPm->GetDoubleParameter(GetFullParmName("Mitochondria/c"), "Length" );}
        
        G4String subComponentName2 = "Mitochondria";
        G4Ellipsoid* gMito = new G4Ellipsoid("gMito", EllA, EllB, EllC);
        G4LogicalVolume* lMito = CreateLogicalVolume(subComponentName2, gMito);
        
        //Randomly distribute mitochondria throughout cell volume
        for (int j = 0; j < NbOfMito; j++){
            auto OverlapCheck = true;
            while (OverlapCheck)
            {
                G4double u = G4UniformRand()*2*pi;
                G4double v = std::acos(2*G4UniformRand()-1);
                G4double dr = G4UniformRand()*(NuclRadius);
                G4double phi = G4UniformRand()*2*pi;
                G4double psi = G4UniformRand()*2*pi;
                G4double x = 0.0;
                G4double y = 0.0;
                G4double z = 0.0;
                
                x = (dr)* std::cos(u) * std::sin(v);
                y = (dr)* std::sin(u) * std::sin(v);
                z = (dr)* std::cos(v);
                
                G4ThreeVector* position = new G4ThreeVector(x,y,z);
                
                G4RotationMatrix* rotm = new G4RotationMatrix();
                
                rotm->rotateX(psi);
                rotm->rotateY(phi);
                
                G4VPhysicalVolume* pMito = CreatePhysicalVolume(subComponentName2, j, true, lMito, rotm, position, pNucleus);
                
                OverlapCheck = pMito->CheckOverlaps();
                
                if (OverlapCheck == true) {
                    pNucleus->GetLogicalVolume()->RemoveDaughter(pMito);
                    pMito = NULL;
                    G4cout << "**** Finding new position for volume " << subComponentName2 << ":" << j <<  " ****" << G4endl;
                }
            }
        }
    } else {
            G4cerr << "Mitochondria Not Defined" << G4endl;
            exit(2);
        }
    }

    InstantiateChildren(fEnvelopePhys);
	
	return fEnvelopePhys;
}
