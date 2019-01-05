/***************************************************************************
 * 
 *           Module :  main.cpp
 *          Program :  cphdShell
 *       Created by :  Darren Muff on 09/12/12016
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *   Program to create a CPHD file from user input. Clues are provided along
 *   the way to ensure that the parameters are consistent within the CPHD file.
 *   The cphd file is a shell of an actual cphd file meaning that it contains no
 *   wideband information (just zeroes sop that it remains the correct size).
 * 
 *   (c) Crown Copyright 2018 Defence Science and Technology Laboratory
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software")
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 * 
 ***************************************************************************/

#include <iostream>
#include <sarclib/sarclib.h>
#include "cphdShell_version.h"
#include "colourCodes.h"
#include "ecef2SceneCoords.hpp"
#include <algorithm>


SPVector ptOnStraightLineTrajectory(SPVector SACentre, SPVector velocity, double tRelToSACentre) ;

typedef enum { SPOTLIGHT, STRIPMAP } SARMode_type;

void banner (){
    printf(" \n");
    printf(DARK GREEN "           cphdShell - A Programme to Create Blank CPHD Shells\n" NORMAL);
    printf(BLUE "                         Version :" RED" %s \n", SHORT_VERSION);
    printf(BLUE "                          Commit :" RED" %s \n", COMMIT);
    printf(BLUE "                Revision: " RED"%s, %s \n",REVISION, VERSION_DATE);
    printf(BLUE "           Copyright (c) 2016 " WHITE"[" BLUE"Dstl" WHITE"]" BLUE". All rights reserved.\n" RESETCOLOR);
    printf(" \n");
    
    return ;
}


int main(int argc, const char * argv[]) {
    
    banner() ;
    
    SPVector Zhat ;    VECT_CREATE(0, 0, 1, Zhat) ;
    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"cphdShell", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status);
    
    
    printf("\n");
    printf("                             GEOMETRY \n");
    printf("                             -------- \n");
    
    SPVector SRPlla ;
    VECT_CREATE(0, 0, 0, SRPlla);
    SRPlla = input_vect("Where is the scene centre in lat/lon/alt?", "SRPlla",
                        "Enter the location on the Earth for the stabilisation reference point \n\
                        (SRP) in decimal degrees and metres", SRPlla) ;
    
    
    double TxRslant = 1000.0;
    TxRslant = input_dbl("What is the slant range at closest approach?", "TxRslant",
                         "This is the range from the transmitter location at the point of\n\
                         closest approach to the scene centre (in metres)", TxRslant);
    
    double TxGraz = 45;
    TxGraz = input_dbl("What is the grazing angle?", "TxGraz",
                       "This is the angle in degrees between the horizontal tangent to the \n\
                       Earth's surface and the direction to the transmitter \n\
                       subtended at the scene centre", TxGraz) ;
    
    double TxAzim = 90 ;
    TxAzim = input_dbl("What is the azimuth angle?", "TxAzim",
                       "This is the angle in degrees subtended at the scene centre in the horizontal \n\
                       plane between the direction of north and the direction to the transmitter.",
                       TxAzim) ;
    
    int TXlookDir = 0 ;
    char *lookstr ;
    lookstr = input_string("What is the look direction (l/r) ?", "TXlookDir",
                          "Which side of the the transmitter platform is the radar \n\
                          transmitting relative to the velocity vector of the platform.\n\
                          (can be 'l','r','left','right')", "left") ;

    if(!(strcmp(lookstr, "left")) || !(strcmp(lookstr,"l") || !(strcmp(lookstr,"L")) )){
        TXlookDir = 0 ;
        printf("Setting look direction to be 'left'\n");
    }else{
        TXlookDir = 1 ;
        printf("Setting look direction to be 'right'\n");
    }
    
    double TXSquint = 0.0 ;
    TXSquint = input_dbl("What is the Transmitter squint angle?", "TXSquint",
                         "This is the angle in degrees from broadside of the transmitting beam.\n\
                         Positive angles are forward from broadside and negative point towards\n\
                         the rear (relative to the velocity vector).", TXSquint) ;
    
    bool modenotset = true;
    SARMode_type sarmode = SPOTLIGHT;
    do{
        char *mode;
        mode = input_string("Sensor mode?", "SARMODE",
                            "Options are \"STRIPMAP\" or \"SPOTLIGHT\"", "SPOTLIGHT");
        std::string SARModeStr  = std::string(mode) ;
        std::transform(SARModeStr.begin(), SARModeStr.end(), SARModeStr.begin(), ::toupper);
        if (SARModeStr == "STRIPMAP") {
            sarmode = STRIPMAP ;
            modenotset = false ;
        }else if (SARModeStr == "SPOTLIGHT"){
            sarmode = SPOTLIGHT ;
            modenotset = false ;
        }else{
            printf("Sensor mode must be either \"STRIPMAP\" or \"SPOTLIGHT\"\n");
        }
    }while (modenotset) ;
    
    bool bistatic = false ;
    bistatic = input_yesno("Is this a bistatic collection?", "bistatic",
                           "If this is yes then there will be a seperate set of questions for \n\
                           the receiving geometry.", bistatic) ;
    if (bistatic) {
        printf("ERROR: Bistatic not yet implemented. Check back in 30 years\n");
        exit(1);
    }
    
    SPVector SAcent ;
    SPVector RvectG, AvectG, RvectS, TxVdir,SRPDir ;
    SPVector SRP ;      latlon2ecef(&SRPlla, &SRP, &status) ;
    SPVector localE ;   VECT_CROSS(Zhat, SRP, localE) ;     VECT_NORM(localE, localE) ;
    SPVector localN ;   VECT_CROSS(SRP, localE, localN) ;   VECT_NORM(localN, localN) ;
    VECT_NORM(SRP, SRPDir) ;
    vectRotateAxis(localN, SRPDir, DEG2RAD(-TxAzim), &RvectG) ;       // Creates the Range vector along the ground
    VECT_CROSS(SRP, RvectG, AvectG) ;   VECT_NORM(AvectG, AvectG) ;     // Creates Azimuth vector along the ground
    vectRotateAxis(RvectG, AvectG, DEG2RAD(-TxGraz), &RvectS) ;       // Creates Range vector in slant plane
    VECT_NORM(RvectS, RvectS); VECT_SCMULT(RvectS, TxRslant, RvectS) ;  // RvectS is now the vector from SRP to SA centre
    VECT_ADD(SRP, RvectS, SAcent) ;                                     // SAcent is now the ECEF coords for the centre of the SA
    SPVector ptRange  ; VECT_SUB(SRP, SAcent, ptRange) ;
    SPVector rVecDir  ; VECT_NORM(ptRange, rVecDir) ;
    SPVector rCross   ; VECT_CROSS(rVecDir, SRPDir, rCross) ; VECT_NORM(rCross, rCross) ;
    SPVector rGrnd    ; VECT_CROSS(SRPDir, rCross, rGrnd) ;
    SPVector rGrndDir ; VECT_NORM(rGrnd, rGrndDir) ;
    
    if(TXlookDir == 0){
        // Left looking
        //
        SPVector tmp ;
        vectRotateAxis(rGrndDir, SRPDir, -((SIPC_pi/2)-DEG2RAD(TXSquint)), &tmp) ;
        VECT_NORM(tmp, TxVdir);
        
        
    }else{
        // Right looking
        SPVector tmp ;
        vectRotateAxis(rGrndDir, SRPDir, ((SIPC_pi/2)-DEG2RAD(TXSquint)), &tmp) ;
        VECT_NORM(tmp, TxVdir);
    }
    
    SPVector SAcent_lla_rb, SAcent_loc_rb, TxVdir_tmp, TxVdir_loc_rb ;
    VECT_ADD(SRP, TxVdir, TxVdir_tmp) ;
    SPVector *points = (SPVector *)sp_malloc(sizeof(SPVector) * 2) ;
    points[0] = SAcent; points[1] = TxVdir_tmp ;
    ecef2SceneCoords(2, points, SRP) ;
    SAcent_loc_rb = points[0]; TxVdir_loc_rb = points[1] ;
    ecef2latlon(&SAcent, &SAcent_lla_rb, &status) ;
    
    printf("                                                Geometry Summary\n");
    printf("              ---------------------------------------------------------------------------------------\n");
    printf("              | lat(deg) lon(deg) alt(m)   |               ECEF                |    Local(Rel SRP)    \n");
    printf("SRP           |     %05.2f,%06.2f,%6.1f    |    %09.1f,%09.1f,%09.1f  |    %4.2f,%4.2f,%4.2f\n",
           SRPlla.x,SRPlla.y,SRPlla.z, SRP.x,SRP.y,SRP.z,0.,0.,0.);
    printf("Syn Ap Cent   |     %05.2f,%06.2f,%08.1f  |    %09.1f,%09.1f,%09.1f  |    %4.2f,%4.2f,%4.2f\n",
           SAcent_lla_rb.x,SAcent_lla_rb.y,SAcent_lla_rb.z, SAcent.x,SAcent.y,SAcent.z, SAcent_loc_rb.x,SAcent_loc_rb.y,SAcent_loc_rb.z);
    printf("Velocity dirn |     -----, ------,------   |         %04.1f,     %04.1f,     %04.1f  |    %4.2f,%4.2f,%4.2f\n",
           TxVdir.x,TxVdir.y,TxVdir.z,TxVdir_loc_rb.x,TxVdir_loc_rb.y,TxVdir_loc_rb.z);
    
    
    printf("\n\n");
    printf("                                               The Transmitter \n");
    printf("                                               --------------- \n");
    
    double fcent = 10e9;
    fcent = input_dbl("What is the centre frequency for the transmitter?", "fcent",
                      "This is the transmitter centre frequency in Hz", fcent) ;
    
    double TxAntLen = 1.0;
    TxAntLen = input_dbl("What is the transmit antenna length?", "TxAntLen",
                         "This is the length in metres in the horizontal plane of the\n\
                         transmit antenna", TxAntLen) ;
    
    double TxAntHgt = 1.0;
    TxAntHgt = input_dbl("What is the transmit antenna height?", "TxAntHgt",
                         "This is the height in metres of the transmit antenna in the direction\n\
                         perpendicular to the range direction and the antenna length direction.",
                         TxAntHgt) ;
    
    double TxBeamHgtAng, TxBeamLenAng, TxBeamHgtDis, TxBeamLenDis, TxBeamHgtDisGrd, wavelength ;
    wavelength      = SIPC_c / fcent ;
    TxBeamHgtAng    = wavelength / TxAntHgt ;
    TxBeamHgtDis    = TxRslant * TxBeamHgtAng ;
    TxBeamHgtDisGrd = TxBeamHgtDis / sin(DEG2RAD(TxGraz)) ;
    TxBeamLenAng    = wavelength / TxAntLen ;
    TxBeamLenDis    = TxRslant * TxBeamLenAng ;
    
    printf("Wavelength : %f m\n", wavelength) ;
    printf("Transmitter beam: azim: %f deg x elev: %f deg \n", DEG2RAD(TxBeamLenAng), DEG2RAD(TxBeamHgtAng));
    printf("Transmitter beam at scene centre: azim: %f m x elev %f m (ground %f m)\n", TxBeamLenDis, TxBeamHgtDis, TxBeamHgtDisGrd) ;
    
    printf("\n");
    printf("                     Pulse Timing \n");
    printf("                     ------------ \n");
    
    double TxVel = 100.0 ;
    TxVel = input_dbl("What is the transmitter platform velocity?", "TxVel",
                      "This is the velocity in m/s of the transmitting platform", TxVel) ;
    SPVector TxVelVec;
    VECT_SCMULT(TxVdir, TxVel, TxVelVec) ;
    printf("Velocity vector in ECEF is %f,%f,%f\n", TxVelVec.x,TxVelVec.y,TxVelVec.z) ;
    // find max and min velocity of beam
    //
    double vmax,vmin,maxAng,minAng,fmin,fmax,Bdopp;
    maxAng = DEG2RAD(TXSquint) + TxBeamLenAng / 2 ;
    minAng = DEG2RAD(TXSquint) - TxBeamLenAng / 2 ;
    vmax   = TxVel * cos((SIPC_pi/2) - maxAng) ;
    vmin   = TxVel * cos((SIPC_pi/2) - minAng) ;
    fmax   = ((SIPC_c + vmax) / SIPC_c) * fcent ;
    fmin   = ((SIPC_c + vmin) / SIPC_c) * fcent ;
    Bdopp  = fmax - fmin ;
    
    printf(" Beam velocity : %f - %f m/s\n",vmin,vmax);
    printf(" Beam frequency : %f - %f Hz\n",fmin,fmax);
    printf(" Doppler bandwidth : %f Hz\n",Bdopp);
    double prfMaxBeam, prfMinBeam;
    prfMinBeam = Bdopp ;
    prfMaxBeam = SIPC_c /(2 * TxBeamHgtDis / tan(DEG2RAD(TxGraz))) ;
    printf(" Minimum PRF to unambiguously sample beam (azim doppler)   : %f Hz\n", prfMinBeam) ;
    printf(" Maximum PRF to unambiguously sample beam (rnge ambiguity) : %f Hz\n", prfMaxBeam) ;
    
    double PRF = 2000;
    double unambigBeamAz ;
    double imageX = 100.0;
    double imageY = 100.0;
    double imageYSlant ;
    
    do {
        PRF = input_dbl("What PRF should be used ?", "PRF",
                        "This is the pulse repetition frequency to use in Hertz", PRF);
        
        unambigBeamAz = ((PRF-(prfMinBeam-PRF))/prfMinBeam) * TxBeamLenDis ;
        if (PRF < prfMinBeam) {
            printf("Unambiguous beam size in azimuth is %f m (PRF limited)\n",unambigBeamAz) ;
        }else{
            printf("Unambiguous beam size in azimuth is %f m (Beam limited)\n",TxBeamLenDis) ;
        }
        
        imageX = input_dbl("How large is the imaged area in Azimuth?", "imageX",
                           "This is the expected length of the spotlight image in the along track direction \n\
                           and is used to guide the ambiguities for the PRF", imageX);
        
        imageY = input_dbl("How large is the imaged area on the ground in Range?", "imageY",
                           "This is the expected size of the spotlight image in the range direction \n\
                           and is used to determine the receive window size", imageY);
        
        if (unambigBeamAz < 0 || ((imageX > unambigBeamAz) & (sarmode == SPOTLIGHT)) ) {
            printf("Warning desired image size in azmuth is not sampled high enough by the PRF to be unambiguous\n");
        }
        
    } while (unambigBeamAz < 0 || ((imageX > unambigBeamAz) & (sarmode == SPOTLIGHT)) );
    
    imageYSlant = imageY * cos(TxGraz) ;
    
    printf("\n");
    printf("                      Resolution  \n");
    printf("                      ----------- \n");
    
    double azRes=1.0,raRes=1.0;
    azRes = input_dbl("What is the required azimuth resolution?", "azRes",
                      "This is the azimuth resolution that the CPHDfile will be built for",azRes) ;
    raRes = input_dbl("What is the required range resolution?", "raRes",
                      "This is the range resolution that the CPHDFile will be built for", raRes) ;
    int groundres = 1;
    groundres = input_yesno("Is this range resolution on the ground?","groundres",
                            "if yes then slant range resolution will be calculated to meet the \n\
                            desired ground resolution using the grazing angle specified", groundres) ;
    
    if (groundres) {
        raRes = raRes * cos(DEG2RAD(TxGraz)) ;
    }
    
    double oversampAz = 1.2, oversampRa=1.2 ;
    oversampAz = input_dbl("What is the Azimuth oversampling factor?","oversampAz",
                           "This is the factor by which the number of pulses required to \n\
                           form the synthetic aperture is multiplied by", oversampAz) ;
    oversampRa = input_dbl("What is the range oversampling factor?", "oversampRa",
                           "This is the factor by which the number of required ADC samples is multiplied by. The CPHD Spec recommends that this should be at least 1.1", oversampRa) ;
    
    double SAduration, SAlen, bandwidth;
    int Npulses;
    
    // SAlen = (TxRslant + (imageYSlant/2)) * wavelength / azRes;
    SAlen = TxRslant * wavelength * oversampAz / (2.0 * azRes * cos(TXSquint)) ;
    SAduration = SAlen / TxVel ;
    if (sarmode == SPOTLIGHT) {
        Npulses = ceil(oversampAz * SAduration * PRF) ;
    }else if (sarmode == STRIPMAP){
        Npulses = oversampRa * (imageX + SAduration * PRF) ;
    }else{
        printf("Error : unknown SAR mode used \n");
        exit(1);
    }
    
    bandwidth = SIPC_c/(2*raRes);
    printf(" Azimuth resolution: %f m\n",azRes);
    printf(" Synthetic aperture duration: %f secs\n", SAduration);
    printf(" Synthetic aperture length : %f m\n", SAlen);
    printf(" Number of pulses in CPHD (inc oversampling) : %d pulses\n",Npulses);
    printf(" Slant range resolution : %f m\n",raRes);
    printf(" Transmit bandwidth is %f MHz\n", bandwidth/1e6 );
    printf(" Pulse transmit start-stop : %f - %f MHz\n", (fcent-(bandwidth/2))/1e6, (fcent+(bandwidth/2))/1e6);
    
    
    printf("\n");
    printf("                      Pulse Sampling  \n");
    printf("                      -------------- \n");
    
    
    printf("The interpulse period for the specified PRF is: %f ms\n", (1/PRF)*1000.0);
    char *pulse_cstr;
    pulse_cstr = (char *)sp_malloc(sizeof(char)*255);
    sprintf(pulse_cstr,"10%%") ;
    pulse_cstr = input_string("What is the desired pulse length (in ms) (or duty cycle )?", "pulselen",
                              "This is the pulse length in milliseconds. A duty cycle can be used if\n\
                              a '%' symbol is placed after the input number", pulse_cstr) ;
    
    char *tok;
    double dutyCycle;
    double pulseLen;
    int len = (int)strlen(pulse_cstr);
    tok = strtok(pulse_cstr, "%");
    if (strlen(tok) != len) {
        dutyCycle = atof(tok) / 100. ;
        pulseLen = dutyCycle / PRF ;
    }else{
        pulseLen  = atof(pulse_cstr) / 1000.0;
        dutyCycle = pulseLen * PRF ;
    }
    printf("Pulse length is %f ms\n",pulseLen*1000) ;
    if (SIPC_c * pulseLen / 2  > TxRslant) {
        printf("Warning Pulse length is larger than slant range for Transmitter - Just saying...\n");
    }
    printf("Duty cycle is %f %%\n",dutyCycle*100);
    printf("Pulse length in flight : %6.3f km\n", SIPC_c*pulseLen/1000);
    
    double gamma = bandwidth/pulseLen;
    printf("Chirp Rate (Gamma) : %e Hz/Sec\n", gamma) ;
    
    double pulseFlightTime = TxRslant * 2 / SIPC_c ;
    double pulsesInFlight = pulseFlightTime * PRF ;
    int nPulseInFlight = floor(pulsesInFlight) ;
    printf("Number of pulses in flight (pulse ambiguity number): %d\n",nPulseInFlight);
    
    double Bif = oversampRa * 2 * gamma * imageY / SIPC_c ;
    printf("IF bandwidth of swath is %f MHz\n", Bif/1e6) ;
    
    double ADCRate = 60; // Hz
    ADCRate = input_dbl("What is the ADC rate in MHz?", "ADCRate",
              "The is the clock speed of the analogue-to-digital converter", ADCRate) ;
    
    ADCRate = ADCRate * 1e6 ;
    if (Bif > ADCRate) {
        printf("Error : ADC rate is not high enough to sample requested swath\n");
        exit(0);
    }
    
    double nSamples ;
    double sampSize = SIPC_c / (2 * ADCRate) ;
    double RxWindow = pulseLen + (2 * imageY / SIPC_c) ;
    nSamples = RxWindow * ADCRate ;
    
    printf("ADC Sample size : %f m\n", sampSize) ;
    printf("Receive window duration : %f ms\n", RxWindow*1000);
    printf("number of samples : %f\n", nSamples);
    
    CPHDHeader hdr ;
    
    time_t rawtime ;
    struct tm * timeinfo ;
    hdr.dateTime = (char *)sp_malloc(sizeof(char) * 80);
    time(&rawtime) ;
    timeinfo = localtime(&rawtime) ;
    strftime(hdr.dateTime,80,"%Y%m%d%H%M%S",timeinfo);
    
    hdr.dateTime = input_string("What DateTime string should be used?", "dateTime",
                 "This is string to use for the dateTime in the CPHD File", hdr.dateTime) ;
    
    hdr.dataSetID = (char *)sp_malloc(sizeof(char) * 255);
    sprintf(hdr.dataSetID, "%s : %05.2f x %05.2f",hdr.dateTime, azRes, raRes);
    
    hdr.dataSetID = input_string("Enter the Dataset ID to use", "datasetID",
                                 "Enter a string to use in the CPHD header file as a dataset ID", hdr.dataSetID);
    
    
    printf("\n");
    printf("                           Output  \n");
    printf("                      -------------- \n");
    char *outCPHDFile;
    outCPHDFile = input_string("Enter name of cphdfile to create", "outCPHDFilename",
                               "The fully qualified pathname of the CPHD File to be crated", "outputCPHDFile.cph") ;
    char *fType;
    bool pass;
    do{
        fType = input_string("What type of CPHD file to create [\'3\' or \'x\']", "CPHDTYPE",
                             "cphd3 is simplest (from SNL), \'x\' is cphdx v0.3", "3") ;
        if(!strcmp(fType, "x")){
            hdr.version = 'x' ;
            pass = true;
        }else if(!strcmp(fType, "3")){
            hdr.version = '3' ;
            pass = true ;
        }else{
            printf("Error : Unsupported CPHD version \'%s\'. should be either \'3\' or \'x\'\n",fType);
            pass = false;
        }
    }while (!pass) ;
    hdr.polarisation = (char **)sp_malloc(sizeof(char **) * 1);
    hdr.polarisation[0] = (char *)sp_malloc(sizeof(char) * 8) ;
    sprintf(hdr.polarisation[0], "VV") ;
    int valid = 1;
    do{
        hdr.polarisation[0] = input_string("Enter polarisation of collection (VV,VH,HV,HH)", "polarisation",
                                           "Polarisation as a string. Must be one of the four options", hdr.polarisation[0]) ;
        if ( (strcmp(hdr.polarisation[0], "VV")) &&  (strcmp(hdr.polarisation[0], "VH")) && (strcmp(hdr.polarisation[0], "HV")) && (strcmp(hdr.polarisation[0], "HH")) ) {
            printf("Error requested polaristion is not one of VV, VH, HV or HH\n");
            valid = 0;
        }else{
            valid = 1;
        }
    }while (valid == 0) ;
    hdr.sensor = (char *)sp_malloc(sizeof(char) * 80) ;
    sprintf(hdr.sensor, "Undefined");
    hdr.sensor = input_string("What sensor name should be used in the cphd file?", "sensorName",
                              "Enter the sensor name to use", hdr.sensor) ;

    printf("Calculating Narrowaband data...\n");
    
    hdr.num_azi = Npulses ;
    hdr.pulses  = (CPHDPulse *)sp_calloc(hdr.num_azi, sizeof(CPHDPulse));
    hdr.antenna_width_az = TxBeamLenAng ;
    hdr.antenna_width_el = TxBeamHgtAng ;
    hdr.chirp_gamma = gamma ;
    hdr.classification = (char *)sp_malloc(sizeof(char) * 10);
    sprintf(hdr.classification,"OFFICIAL") ;
    hdr.clock_speed = ADCRate ;
    hdr.data_type = ITYPE_CMPL_FLOAT ;
    hdr.deskewed  = 0;
    hdr.fixedSRP  = 1;
    hdr.freq_centre = fcent ;
    hdr.geometry = (char *)sp_malloc(sizeof(char)*80);
    sprintf(hdr.geometry,"Monostatic");
    hdr.interleaved = 0;
    hdr.mode = (char *)sp_malloc(sizeof(char)*80);
    if (sarmode == SPOTLIGHT) {
        sprintf(hdr.mode,"SPOTLIGHT");
    }else if (sarmode == STRIPMAP){
        sprintf(hdr.mode,"STRIPMAP");
    }
    hdr.nchan = 1;
    hdr.num_azi = Npulses ;
    hdr.pulse_length = pulseLen ;

    if (hdr.version == 'x') {
        // CPHDX v0.3 is ambiguous in terms of pulse length and record window
        // so the pulse length is commonly set to 0.95 the record window
        //
        hdr.TOASaved = hdr.pulse_length / 0.95 ;
        hdr.nsamp = hdr.TOASaved * hdr.clock_speed ;

    }else{
        // CPHD3 allows us to store the clock speed and pulse length
        //
        hdr.TOASaved = hdr.nsamp / hdr.clock_speed;
        hdr.nsamp = nSamples ;
    }
    hdr.phaseSgn = -1 ;
    hdr.byte_swap = 0 ;
    hdr.num_nb_items = 10 ;
    hdr.pos_nb[0] = CPHD_NB_ChannelNumber ;
    hdr.pos_nb[1] = CPHD_NB_VectorNumber  ;
    hdr.pos_nb[2] = CPHD_NB_SRP           ;
    hdr.pos_nb[3] = CPHD_NB_TxPos         ;
    hdr.pos_nb[4] = CPHD_NB_RcvPos        ;
    hdr.pos_nb[5] = CPHD_NB_TxTime        ;
    hdr.pos_nb[6] = CPHD_NB_RcvTime       ;
    hdr.pos_nb[7] = CPHD_NB_Fx0           ;
    hdr.pos_nb[8] = CPHD_NB_FxStepSize    ;
    hdr.pos_nb[9] = CPHD_NB_AmpSF0        ;
    
    hdr.num_wb_items = 3 ;
    hdr.pos_wb[0] = CPHD_WB_ChannelNumber ;
    hdr.pos_wb[1] = CPHD_WB_VectorNumber  ;
    hdr.pos_wb[2] = CPHD_WB_Samples       ;
    
    for (int i =0; i<Npulses; ++i) {
        // pulse number
        //
        hdr.pulses[i].pulse_number = i ;
        
        // TxTime
        //
        double TxTime = i / PRF ;
        hdr.pulses[i].sat_tx_time  = TxTime ;
        
        // SRP
        //
        SPVector pulseSRP;
        if(sarmode == STRIPMAP){
            pulseSRP = ptOnStraightLineTrajectory(SRP, TxVelVec, TxTime - ((Npulses/2) / PRF)) ;
        }else{
            pulseSRP = SRP ;
        }
        hdr.pulses[i].srp = pulseSRP ;
        
        // TxPos
        //
        hdr.pulses[i].sat_ps_tx = ptOnStraightLineTrajectory(SAcent, TxVelVec, TxTime - ((Npulses/2) / PRF)) ;
        
        // RxTime - Assume Start/Stop Approx. and same as Tx
        //
        hdr.pulses[i].sat_rx_time = TxTime ;
        
        // RxPos - Assume Start/Stop Approx. and same as Tx
        //
        hdr.pulses[i].sat_ps_rx = hdr.pulses[i].sat_ps_tx;
        
        // fx_stap_size
        //
        hdr.pulses[i].fx_step_size = hdr.chirp_gamma / hdr.clock_speed ;
        
        // fx0
        //
        hdr.pulses[i].fx0 = hdr.freq_centre - ((hdr.nsamp/2) * hdr.pulses[i].fx_step_size) ;
        
        // amp_sf0
        //
        hdr.pulses[i].amp_sf0 = 1 ;
        
        // fx1
        //
        hdr.pulses[i].fx1 = hdr.freq_centre - (hdr.chirp_gamma * pulseLen / 2) ;
        
        // fx2
        //
        hdr.pulses[i].fx2 = hdr.freq_centre + (hdr.chirp_gamma * pulseLen / 2) ;
        
        // Wideband data - empty buffer
        //
        hdr.pulses[i].data.cmpl_f = (SPCmplx *)sp_malloc(sizeof(SPCmplx) * hdr.nsamp) ;
        for (int j=0; j<hdr.nsamp; ++j) {
            hdr.pulses[i].data.cmpl_f[j].r = 0.0f ;
            hdr.pulses[i].data.cmpl_f[j].i = 0.0f ;
        }
    }
    
    hdr.fp = fopen(outCPHDFile, "w") ;
    if (hdr.fp == NULL) {
        fprintf(stderr, "Failed to open file %s\n", outCPHDFile);
        perror("Error opening file");
        exit(1);
    }
    printf("Writing CPHD File \"%s\"....",outCPHDFile);
    writeCPHDFile(&hdr, &status) ;
    CHECK_STATUS(&status) ;
    fclose(hdr.fp);
    printf("...Done\n");
    
    im_close_lib(&status);
    
    return 0;
}

SPVector ptOnStraightLineTrajectory(SPVector SACentre, SPVector velocity, double tRelToSACentre){
    SPVector pos, ans;
    VECT_SCMULT(velocity, tRelToSACentre, pos) ;
    VECT_ADD(SACentre, pos, ans) ;
    return (ans);
}


