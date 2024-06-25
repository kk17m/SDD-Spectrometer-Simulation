# SDD-Spectrometer-Simulation
A Monte Carlo simulation code of silicon drift detectors (SDDs) for X-ray fluorescence and spectrometry applications. The simulation includes energy-dispersive photon counting detectors intended for spectrometry and molecular imaging. The detector is based on the Amptek XR100 FastSDD X-ray  detector (Amptek Inc., Bedford, MA, USA).

-------------------------------------------------------------------

     =========================================================
        Silicon drift detector (SDD) simulation on Geant4
     =========================================================

                                SDD
                            -----------
-------------------------------------------------------------------

 This simulation code implements a silicon drift detector (SDD) for x-ray fluorescence and
 spectrometry applications. The detector parameters are based on Amptek XR100 FastSDD X-ray
 detector (Amptek Inc., Bedford, MA, USA). 

 Author: Kunal Kumar (kunal.kumar@ovgu.de)

 Citation: Kumar K, Fachet M, Hoeschen C. High-Spatial-Resolution Benchtop X-ray Fluorescence Imaging
 through Bragg-Diffraction-Based Focusing with Bent Mosaic Graphite Crystals: A Simulation Study.
 Int J Mol Sci. 2024 Apr 26;25(9):4733. doi: 10.3390/ijms25094733. PMID: 38731956; PMCID: PMC11083219.

 DOI: [10.3390/ijms25094733](https://doi.org/10.3390/ijms25094733); PMID: [38731956](https://pubmed.ncbi.nlm.nih.gov/38731956/); PMCID: [PMC11083219](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11083219/)

 -------------------------------------------------------------------
	
 1- GEOMETRY DEFINITION
	
   The geometry is constructed in the SDDDetectorConstruction class. The parameters for geometry are taken
   from the user manual (Amptek Silicon Drift Detector XR-100FastSDD). [Available online](https://www.amptek.com/-/media/ametekamptek/documents/resources/products/user-manuals/xr100-1mm-fastsdd-user-manual-b4.pdf?la=en&revision=24e8eb09-6164-48ba-8336-e572f84bf5c1)
   		
 2- PHYSICS LIST
   
   The physics list is based on the Geant4 TestEm5 example (electromagnetic/TestEm5/src/PhysicsList.cc).
   Low-energy EM physics is implemented using the Penelope model.
  	 
 3- PRIMARY GENERATOR
  
   The primary generator is defined in the SDDPrimaryGeneratorAction class.
   The default kinematics is a 45 keV monoenergetic gamma pencil beam. See
   SDD.in.
     
 4- SDD DETECTOR MODEL

   The detector model is based on analytical equations implementing Fano and electronic noise. The parameters
   are obtained from the Amptek user manual and application Notes (Amptek Silicon Drift Detector XR-100FastSDD). [Available online](https://www.amptek.com/-/media/ametekamptek/documents/resources/application-notes/high-sensitivity-detectors-for-xrf.pdf?la=en&revision=9d04dd37-c2ea-4f89-ad58-55579a8574b1)
   
 5- HOW TO RUN

    - Execute SDD in the 'interactive mode' with visualization:
        % ./SDD
      and type in the commands from SDD.in line by line:  
        Idle> /control/verbose 2
        Idle> /tracking/verbose 1
        Idle> /run/beamOn 10 
        Idle> ...
        Idle> exit
      or
        Idle> /control/execute SDD.in
        ....
        Idle> exit

    - Execute SDD in the 'batch' mode from macro files 
      (without visualization)
        % ./SDD SDD.in
        % ./SDD SDD.in > SDD.out

	
