/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/* $Id$ */

//_________________________________________________________________________
// Implementation version v0 of EMCAL Manager class 
// An object of this class does not produce hits nor digits
// It is the one to use if you do not want to produce outputs in TREEH or TREED
// This class places a Geometry of the EMCAL in the ALICE Detector as defined in AliEMCALGeometry.cxx                 
//*-- Author: Yves Schutz (SUBATECH)
//*-- and   : Sahal Yacoob (LBL / UCT)

// This Version of AliEMCALv0 reduces the number of volumes placed in XEN1 (the envelope) to less than five hundred
// The Envelope is Placed in Alice, And the Aluminium layer. Mini envelopes (XU) are then placed in XEN1.
// Each mini envelope contains 2 scintillator, and 2 lead layers, except the last one which contains just one scintillator layer.
// At the moment I cannot place the 36 and above layers in the mini envelopes so all layers are still placed in XEN1


// --- ROOT system ---

#include "TPGON.h"
#include "TTUBS.h"
#include "TNode.h"
#include "TRandom.h"
#include "TGeometry.h"

// --- Standard library ---

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Rstrstream.h>
#include <Riostream.h>

// --- AliRoot header files ---

#include "AliEMCALv0.h"
#include "AliEMCALGeometry.h"
#include "AliConst.h"
#include "AliRun.h"

ClassImp(AliEMCALv0)

//______________________________________________________________________
AliEMCALv0::AliEMCALv0(const char *name, const char *title):
  AliEMCAL(name,title)
{
  // ctor : title is used to identify the layout
  GetGeometry() ; 

}

//______________________________________________________________________
void AliEMCALv0::BuildGeometry()
{
    // Display Geometry for display.C

    const Int_t kColorArm1   = kBlue ;

    AliEMCALGeometry * geom = GetGeometry() ; 

    // Define the shape of the Calorimeter 
    TNode * top = gAlice->GetGeometry()->GetNode("alice") ;
    new TTUBS("Envelop1", "Tubs that contains arm 1", "void", 
	      geom->GetEnvelop(0),     // rmin 
	      geom->GetEnvelop(1) +30 ,     // rmax
	      geom->GetEnvelop(2)/2.0, // half length in Z
	      geom->GetArm1PhiMin(),   // minimun phi angle
	      geom->GetArm1PhiMax()    // maximun phi angle
	);

    // Place the Node
    top->cd();
    TNode * envelop1node = new TNode("Envelop1", "Arm1 Envelop", "Envelop1"
				     ,0., 0., 0., "") ;
    envelop1node->SetLineColor(kColorArm1) ;
    fNodes->Add(envelop1node) ;
}

//______________________________________________________________________
void AliEMCALv0::CreateGeometry()
{
  // Create the EMCAL geometry for Geant
  // Geometry of a tower
  //|-----------------------------------------------------| XEN1
  //| |                                                 | |
  //| |    Al thickness = GetAlFrontThickness()         | |
  //| |                                                 | |
  //| |                                                 | |
  //| |                                                 | |
  //|  -------------------------------------------------  |
  //| |    Air Gap = GetGap2Active()                    | |
  //| |                                                 | |
  //|  -------------------------------------------------  |
  //| |    XU0 : XPST (PreShower e = GetPreSintThick() )| |
  //|  -------------------------------------------------  |
  //| |    XU0 : XPBX (PreShower e = GetPbRadThick() )  | |
  //|  -------------------------------------------------  |
  //| |    XU0 : XPST (PreShower e = GetPreSintThick() )| |
  //|  -------------------------------------------------  |
  //| |    XU0 : XPBX (PreShower e = GetPbRadThick() )  | |
  //|  -------------------------------------------------  |
  //| |    XU1 : XPST (Tower e = GetFullSintThick() )   | |
  //|  -------------------------------------------------  |
  //| |    XU1 : XPBX (Tower e = GetPbRadThick() )      | |
  //|  -------------------------------------------------  |
  //| |    XU1 : XPST (Tower e = GetFullSintThick()     | |
  //|  -------------------------------------------------  |
  //| |    XU1 : XPBX (Tower e = GetPbRadThick() )      | |
  //|  -------------------------------------------------  |
  //|    etc .....                                        |
  //|  -------------------------------------------------  |
  //| |    XU10 : XPST (Tower e = GetFullSintThick() )  | |
  //|-----------------------------------------------------|
 
    Float_t etamin,etamax;
    Float_t *dum=0;

    AliEMCALGeometry * geom = GetGeometry() ; 

    if(!(geom->IsInitialized())){
	Error("CreateGeometry","EMCAL Geometry class has not been set up.");
    } // end if

    // Get pointer to the array containing media indices
    Int_t *idtmed = fIdtmed->GetArray() - 1599 ;

    Int_t idrotm = 1;
    AliMatrix(idrotm, 90.0, 0., 90.0, 90.0, 0.0, 0.0) ;


    // Create the EMCAL Mother Volume (a polygone) within which to place the Detector and named XEN1 

    Float_t envelopA[10];
    envelopA[0] = geom->GetArm1PhiMin();                         // minimum phi angle
    envelopA[1] = geom->GetArm1PhiMax() - geom->GetArm1PhiMin(); // angular range in phi
    envelopA[2] = geom->GetNPhi();                               // number of sections in phi
    envelopA[3] = 2;                                             // 2 z coordinates
    envelopA[4] = geom->ZFromEtaR(geom->GetEnvelop(0),
				   geom->GetArm1EtaMin());       // z coordinate 1
    envelopA[5] = geom->GetEnvelop(0) ;                          // rmin at z1
    envelopA[6] = geom->GetEnvelop(1) ;                          // rmax at z1
    envelopA[7] = geom->ZFromEtaR(geom->GetEnvelop(0),
				  geom->GetArm1EtaMax());        // z coordinate 2
    envelopA[8] = envelopA[5] ;                                  // radii are the same.
    envelopA[9] = envelopA[6] ;                                  // radii are the same.

    gMC->Gsvolu("XEN1", "PGON ", idtmed[1599], envelopA, 10) ;   // Polygone filled with air 

    // Position the EMCAL Mother Volume (XEN1) in Alice (ALIC)  

    gMC->Gspos("XEN1", 1, "ALIC", 0.0, 0.0, 0.0, idrotm, "ONLY") ;
    
    if (fDebug) 
      Info("CreateGeometry","rXEN1 = %f, %f\n", envelopA[5], envelopA[6]); 

    // Create mini-envelopes which will contain the PreShower scintillator-Lead-scintillator-lead (XU0) 
   
    TString label = "XU0";

    envelopA[5] = envelopA[5] + geom->GetGap2Active() 
      + geom->GetAlFrontThickness();                              // rmin at z1
    envelopA[4] = geom->ZFromEtaR(envelopA[5],
				  geom->GetArm1EtaMin());         // z coordinate 1
    envelopA[7] = geom->ZFromEtaR(envelopA[5],
				  geom->GetArm1EtaMax());         // z coordinate 2
    envelopA[6] = envelopA[5] + 2 * (geom->GetPreSintThick()
      + geom->GetPbRadThick() );                                  // rmax at z1
    envelopA[8] = envelopA[5] ;                                   // radii are the same.
    envelopA[9] = envelopA[6] ;                                   // radii are the same.
    
    gMC->Gsvolu(label.Data(), "PGON", idtmed[1599], envelopA, 10);// Polygone filled with air 

    // Position XU0 in XEN1
  
    gMC->Gspos(label.Data(), 1, "XEN1", 0.0, 0.0, 0.0, idrotm, "ONLY");

    if (fDebug) 
      Info("CreateGeometry","rXU0 = %f, %f\n", envelopA[5], envelopA[6]); 

    // Create mini-envelopes which will contain the Tower scintillator-Lead-scintillator-lead (XU1 -> XU9)
    
    Float_t tseg = geom->GetFullSintThick()+geom->GetPbRadThick(); // thickness of scintillator+Pb
 
    Int_t i ; 
    for (i = 1; i < ((geom->GetNLayers()-2)/2) + 1 ; i++ ){
	label = "XU" ;
	label += i ;
	envelopA[5] = envelopA[6] ;                                   // rmin at z1
	envelopA[4] = geom->ZFromEtaR(envelopA[5],
				      geom->GetArm1EtaMin());         // z coordinate 1
	envelopA[7] = geom->ZFromEtaR(envelopA[5],
				      geom->GetArm1EtaMax());         // z coordinate 2
	envelopA[6] = envelopA[5] + 2 * tseg ;                        // rmax at z1
	envelopA[8] = envelopA[5] ;                                   // radii are the same.
	envelopA[9] = envelopA[6] ;                                   // radii are the same.
 
	gMC->Gsvolu(label.Data(), "PGON", idtmed[1599], envelopA, 10);// Polygone filled with air 

	// Position XUi in XEN1
	
	gMC->Gspos(label.Data(), 1, "XEN1", 0.0, 0.0, 0.0, idrotm, "ONLY") ;

	if (fDebug) 
	  Info("CreateGeometry","rXEN%d = %f, %f\n", i, envelopA[5], envelopA[6]); 

    } // end  i
 
  
    // Create one mini-envelope which will contain the last Tower scintillator (XU(nlayers-1)/2)

    tseg = geom->GetFullSintThick() ;
    label = "XU" ;
    label += i ;
    envelopA[5] = envelopA[6] ;                                   // rmin at z1
    envelopA[4] = geom->ZFromEtaR(envelopA[5],
				  geom->GetArm1EtaMin());         // z coordinate 1
    envelopA[7] = geom->ZFromEtaR(envelopA[5],
				  geom->GetArm1EtaMax());         // z coordinate 2
    envelopA[6] = envelopA[5] + tseg ;                            // rmax at z1
    envelopA[8] = envelopA[5] ;                                   // radii are the same.
    envelopA[9] = envelopA[6] ;                                   // radii are the same.

    gMC->Gsvolu(label.Data(), "PGON", idtmed[1599], envelopA, 10); // Polygone filled with air

    // Position XU10 in XEN1
  
    gMC->Gspos(label.Data(), 1, "XEN1", 0.0, 0.0, 0.0, idrotm, "ONLY") ;
  
    Info("CreateGeometry","rXEN%d = %f, %f\n", i, envelopA[5], envelopA[6]); 
  
    // Create the shapes of active material (LEAD/Aluminium/Scintillator)
    // to be placed
    Float_t envelopB[10]; // First Layer of Aluminium
    Float_t envelopC[10]; // Scintillator Layers
    Float_t envelopD[10]; // Lead Layers

    envelopC[0] = envelopD[0] = envelopB[0] = envelopA[0] ;  // starting position in Phi
    envelopC[1] = envelopD[1] = envelopB[1] = envelopA[1] ;  // angular range in phi          
    envelopC[2] = envelopD[2] = envelopB[2] = envelopA[2] ;  // number of sections in Phi
    envelopD[3] = envelopC[3] = envelopB[3] = envelopA[3] ;  // 2 z coordinates

    envelopB[4] = geom->ZFromEtaR(geom->GetEnvelop(0),
				  geom->GetArm1EtaMin());   // z co-ordinate 1
    envelopB[5] = geom->GetEnvelop(0) ;                     // rmin at z1
    envelopB[6] = envelopB[5] + geom->GetAlFrontThickness();// rmax at z1
    envelopB[7] = geom->ZFromEtaR(geom->GetEnvelop(0),
				  geom->GetArm1EtaMax());   // z co-ordinate 2
    envelopB[8] = envelopB[5] ;                             // radii are the same.
    envelopB[9] = envelopB[6] ;                             // radii are the same.

    // Define active volumes completely
    
    gMC->Gsvolu("XALU", "PGON", idtmed[1602], envelopB, 10); // PGON filled with Al
    
    gMC->Gspos("XALU", 1, "XEN1", 0.0, 0.0, 0.0 , idrotm, "ONLY") ; // Position Aluminium Layer in XEN1

    gMC->Gsvolu("XPST", "PGON", idtmed[1601], dum, 0);      // PGON filled with Scintillator (shape to be defined by GSPOSP)
  
    gMC->Gsvolu("XPBX", "PGON", idtmed[1600], dum, 0);      // PGON filled with Lead (shape to be defined by GSPOSP)
  
    gMC->Gsdvn("XPHI", "XPST", geom->GetNPhi(), 2);         // Divide eta section of scintillators into phi segments.
 
    // Position alternatively scintillator and  Lead Layers in XUi.

    envelopD[6] = envelopB[6] + geom->GetGap2Active() ;// gap between Al layer and XU0

    for (int i = 0; i < geom->GetNLayers() ; i++ ){
	label = "XU" ;
        label += static_cast<Int_t> (i/2)  ; // we will place two layers (i = one layer) in each mini envelope)	
        envelopC[5] = envelopD[6] ; //rmin
	envelopC[6] = envelopD[6] + ((i > 1)  ? geom->GetFullSintThick() : 
				     geom->GetPreSintThick());//rmax larger for first two layers (preshower)
	envelopC[8] = envelopD[6] ; //rmin
	envelopC[9] = envelopD[6] + ((i > 1 ) ? geom->GetFullSintThick() :
				     geom->GetPreSintThick());//rmax larger for first two layers (preshower)
	for (int j =0; j < (geom->GetNEta()) ; j++){
	    etamin = geom->GetArm1EtaMin()+
		(j*geom->GetDeltaEta());
	    etamax = geom->GetArm1EtaMin()+
		((j+1)*geom->GetDeltaEta());
	    envelopC[4] = geom->ZFromEtaR(envelopD[6],etamin); //z begin  
	    envelopC[7] = geom->ZFromEtaR(envelopD[6],etamax);// z end 

	    gMC->Gsposp("XPST",1+j+i*(geom->GetNEta()), label.Data(), 
			0.0, 0.0, 0.0 , idrotm, "ONLY", envelopC, 10); // Position and define layer
	} // end for j

	if (i < (geom->GetNLayers()-1)){ // except for the last XU which contains only one scintillator layer
	    envelopD[5] = envelopC[6] ; //rmin
	    envelopD[6] = envelopC[6] + geom->GetPbRadThick();  //rmax
	    envelopD[8] = envelopC[6] ; //rmin
	    envelopD[9] = envelopC[6] + geom->GetPbRadThick();  //rmax
	    for (int j =0; j < (geom->GetNEta()) ; j++){
		etamin = geom->GetArm1EtaMin()+
		    (j*geom->GetDeltaEta());
		etamax = geom->GetArm1EtaMin()+
		    ((j+1)*geom->GetDeltaEta());
		envelopD[4] = geom->ZFromEtaR(envelopC[6],etamin);//z begin  
		envelopD[7] = geom->ZFromEtaR(envelopC[6],etamax);// z end

		// Position and Define Layer

		gMC->Gsposp("XPBX",1+j+i*(geom->GetNEta()), label.Data(), 
			    0.0, 0.0, 0.0 , idrotm, "ONLY", envelopD, 10);
	    } // end for j
	} // end if i
    }  // for i
}

//______________________________________________________________________
void AliEMCALv0::Init(void)
{
    // Just prints an information message
  
  if(fDebug) { 
    TString message("\n") ; 
    message += "*****************************************" ;
    
    // Here the EMCAL initialisation code (if any!)
    
    AliEMCALGeometry * geom = GetGeometry() ; 
    
    if (geom!=0) {   
      message += "AliEMCAL " ; 
      message += Version() ; 
      message += "EMCAL geometry intialized for " ; 
      message += geom->GetName()  ;
    }
    else {
      message += "AliEMCAL " ; 
      message += Version() ;  
      message += "EMCAL geometry initialization failed !" ; 
    }
    message += "*****************************************" ;
    Info("Init", message.Data() ) ; 
  }
}
