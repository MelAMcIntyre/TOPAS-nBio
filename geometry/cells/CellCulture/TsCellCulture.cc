// Component for TsCellCulture
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

// A simple cell culture consisting of random spherical cells.

#include "TsCellCulture.hh"

#include "TsParameterManager.hh"
#include "G4VPhysicalVolume.hh"

#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"

#include "G4VSolid.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "globals.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <sstream>
#include <vector>
#include <chrono>

#define G4endl std::endl
using namespace std;

TsCellCulture::TsCellCulture(TsParameterManager *pM, TsExtensionManager *eM, TsMaterialManager *mM, TsGeometryManager *gM,
                             TsVGeometryComponent *parentComponent, G4VPhysicalVolume *parentVolume, G4String &name) : TsVGeometryComponent(pM, eM, mM, gM, parentComponent, parentVolume, name)
{
    ResolveParameters();
}

TsCellCulture::~TsCellCulture()
{
    ;
}

void TsCellCulture::ResolveParameters()
{
    // Half length of container in which cells are placed
    HLX = fPm->GetDoubleParameter(GetFullParmName("Container_HLX"), "Length");
    HLY = fPm->GetDoubleParameter(GetFullParmName("Container_HLY"), "Length");
    HLZ = fPm->GetDoubleParameter(GetFullParmName("Container_HLZ"), "Length");

    // Read Radius and Cell Nucleus
    CellRadius = fPm->GetDoubleParameter(GetFullParmName("CellRadius"), "Length");
    NbOfCells = fPm->GetIntegerParameter(GetFullParmName("NumberOfCells"));
    NuclRadius = fPm->GetDoubleParameter(GetFullParmName("NucleusRadius"), "Length");
}

G4VPhysicalVolume *TsCellCulture::Construct()
{
    BeginConstruction();

    //***********************************************************************
    //              Envelope Geometry : Rectangular container
    //***********************************************************************

    G4Box *gBox = new G4Box(fName, HLX, HLY, HLZ);
    fEnvelopeLog = CreateLogicalVolume(gBox);
    fEnvelopePhys = CreatePhysicalVolume(fEnvelopeLog);

    //***********************************************************************
    //              Cell geometry : spherical
    //***********************************************************************
    // Cell geometry
    G4Orb *gCell = new G4Orb("cell", CellRadius);
    G4LogicalVolume *lCell = CreateLogicalVolume(gCell);

    //***********************************************************************
    // Optional : include a nucleus and/or cylinders in the cell
    //***********************************************************************

    // Nucleus
    G4String subComponentName1 = "Nucleus";
    G4Orb *gNucleus = new G4Orb("gNucleus", NuclRadius);
    G4LogicalVolume *lNucleus = CreateLogicalVolume(subComponentName1, gNucleus);

    std::ofstream cellfile("cell_positions.csv");
    std::ofstream dnafile("dna_positions.csv");

    G4bool placeCylinders = fPm->GetBooleanParameter(GetFullParmName("Cylinders/DoRandomPlacement"));

    // Randomly place cells in the volume
    for (int j = 0; j < NbOfCells; j++)
    {

        G4bool cellPlaced = false;
        while (cellPlaced == false)
        {

            G4double phi = 0;
            G4double psi = 0;
            G4double x = 0.0;
            G4double y = 0.0;
            G4double z = 0.0;

            x = (2 * G4UniformRand() - 1) * (HLX - CellRadius);
            y = (2 * G4UniformRand() - 1) * (HLY - CellRadius);
            z = (2 * G4UniformRand() - 1) * (HLZ - CellRadius);

            G4ThreeVector *position = new G4ThreeVector(x, y, z);
            G4ThreeVector *posNucl = new G4ThreeVector(0 * mm, 0 * mm, 0 * mm);

            G4RotationMatrix *rotm = new G4RotationMatrix();

            rotm->rotateX(psi);
            rotm->rotateY(phi);

            G4VPhysicalVolume *pCell = CreatePhysicalVolume("Cell", j, true, lCell, rotm, position, fEnvelopePhys);

            G4bool OverlapCheck = pCell->CheckOverlaps();

            if (OverlapCheck == false)
            {

                G4VPhysicalVolume *pNucleus = CreatePhysicalVolume("Nucleus", j, true, lNucleus, rotm, posNucl, pCell);

                // write cell coordinates to csv
                cellfile << x << "," << y << "," << z << "\n";
                //*******************************
                // Subcomponent: Cylinders
                //*******************************

                G4String name = GetFullParmName("Cylinders/NumberOfCylinders");
                if (fPm->ParameterExists(name))
                {

                    // Number of Cylinders
                    const G4int NbOfCyl = fPm->GetIntegerParameter(GetFullParmName("Cylinders/NumberOfCylinders"));

                    // Radius and Height of cylinder/DNA segment (default values if none are specified)
                    G4double radius = 1.15 * nanometer;
                    G4double height = 6.8 * nanometer;

                    // Read if they are specified
                    name = GetFullParmName("Cylinders/Radius");
                    if (fPm->ParameterExists(name))
                    {
                        radius = fPm->GetDoubleParameter(GetFullParmName("Cylinders/Radius"), "Length");
                    }

                    name = GetFullParmName("Cylinders/Height");
                    if (fPm->ParameterExists(name))
                    {
                        height = fPm->GetDoubleParameter(GetFullParmName("Cylinders/Height"), "Length");
                    }

                    G4String subComponentName2 = "Cylinders";
                    G4Tubs *gCyl = new G4Tubs("gCyl", 0.0 * micrometer, radius, height, 0.0 * deg, 360.0 * deg);
                    G4LogicalVolume *lCyl = CreateLogicalVolume(subComponentName2, gCyl);

                    if (placeCylinders == true)
                    {
                        // Randomly distribute Cylinders throughout cell volume
                        for (int i = 0; i < NbOfCyl; i++)
                        {

                            G4String cylinderPlaced = false;
                            while (cylinderPlaced == false)
                            {
                                G4double u = G4UniformRand() * 2 * pi;
                                G4double v = std::acos(2 * G4UniformRand() - 1);
                                G4double dr = G4UniformRand() * (NuclRadius);
                                G4double phiCyl = G4UniformRand() * 2 * pi;
                                G4double psiCyl = G4UniformRand() * 2 * pi;
                                G4double xCyl = 0.0;
                                G4double yCyl = 0.0;
                                G4double zCyl = 0.0;

                                xCyl = (dr)*std::cos(u) * std::sin(v);
                                yCyl = (dr)*std::sin(u) * std::sin(v);
                                zCyl = (dr)*std::cos(v);

                                G4ThreeVector *positionCyl = new G4ThreeVector(xCyl, yCyl, zCyl);

                                G4RotationMatrix *rotmCyl = new G4RotationMatrix();

                                rotmCyl->rotateX(psiCyl);
                                rotmCyl->rotateY(phiCyl);

                                G4VPhysicalVolume *pCyl = CreatePhysicalVolume(subComponentName2, i, true, lCyl, rotmCyl, positionCyl, pNucleus);

                                G4bool OverlapCheckCyl = pCyl->CheckOverlaps();

                                if (OverlapCheckCyl == true)
                                {
                                    pNucleus->GetLogicalVolume()->RemoveDaughter(pCyl);
                                    G4cout << "**** Finding new position for volume " << subComponentName2 << ":" << i << " ****" << G4endl;
                                    cylinderPlaced = false;
                                }
                                else
                                {
                                    dnafile << j << "," << i << "," << xCyl << "," << yCyl << "," << zCyl << "," << phiCyl << "," << psiCyl << "\n";
                                    cylinderPlaced = true;
                                }
                            }
                        }
                    }
                    else
                    {
                        // Use pre-placed cylinders
                        ifstream in("dna_positions_test1.csv");

                        string line, field;

                        vector<vector<string>> array; // the 2D array
                        vector<string> v;             // array of values for one line only

                        while (getline(in, line)) // get next line in file
                        {
                            v.clear();
                            stringstream ss(line);

                            while (getline(ss, field, ',')) // break line into comma delimitted fields
                            {
                                v.push_back(field); // add each field to the 1D array
                            }
                            array.push_back(v); // add the 1D array to the 2D array
                        }

                        // Randomly distribute Cylinders throughout cell volume
                        for (size_t i = 0; i < array.size(); ++i)
                        {

                            G4double phiCyl = std::stod(array[i][5]);
                            G4double psiCyl = std::stod(array[i][6]);
                            G4double xCyl = std::stod(array[i][2]);
                            G4double yCyl = std::stod(array[i][3]);
                            G4double zCyl = std::stod(array[i][4]);

                            G4ThreeVector *positionCyl = new G4ThreeVector(xCyl, yCyl, zCyl);
                            G4RotationMatrix *rotmCyl = new G4RotationMatrix();

                            rotmCyl->rotateX(psiCyl);
                            rotmCyl->rotateY(phiCyl);

                            G4VPhysicalVolume *pCyl = CreatePhysicalVolume(subComponentName2, i, true, lCyl, rotmCyl, positionCyl, pNucleus);

                            G4cout << "**** Placing " << subComponentName2 << ":" << i << " ****" << G4endl;
                        }
                    }
                }
                cellPlaced = true;
                break;
            }
            else
            {
                delete pCell;
                cellPlaced = false;
                G4cout << "**** Finding new position for volume Cell : " << j << " ****" << G4endl;
            }
        }
    }

    dnafile.close();
    cellfile.close();
    InstantiateChildren(fEnvelopePhys);

    return fEnvelopePhys;
}
