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

#ifndef TsSphericalCell_hh
#define TsSphericalCell_hh

#include "TsVGeometryComponent.hh"


class TsSphericalCell : public TsVGeometryComponent
{    
public:
	TsSphericalCell(TsParameterManager* pM, TsExtensionManager* eM, TsMaterialManager* mM, TsGeometryManager* gM,
				  TsVGeometryComponent* parentComponent, G4VPhysicalVolume* parentVolume, G4String& name);
	~TsSphericalCell();
	
	G4VPhysicalVolume* Construct();
    
    void ResolveParameters();
    
private:
    G4double CellRadius;

};

#endif
