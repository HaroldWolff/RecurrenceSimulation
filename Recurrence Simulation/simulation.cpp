#include "simulation.hpp"

using namespace std;
ofstream Ostream;
Rndm* rndm = new Rndm(time(NULL));

void GenerateOutputfile(int TOT, int i, int j)
{
	/// *** Change path to outputfile here ***

	string outputFileName = "C:/.... path .../simulation_I_" + boost::lexical_cast<string>(i) + "_J_" + boost::lexical_cast<string>(j) + ".csv";
	cout << outputFileName << endl;
    Ostream.open(outputFileName.c_str());

    // add analysis conditions and parameter names
    Ostream << ", All,=AANTAL(C8:C" << 7+TOT << "),,,,,,,, Poly,=AANTAL(L8:L" << 7+TOT << "),,,,,,,, True Oligo,=AANTAL(U8:U" << 7+TOT << "),,,,,,,, FP Oligo,=AANTAL(AC8:AC" << 7+TOT << ")," << endl;
    Ostream << ", R, G, Mtotal, RecurrencesDetected, symptomaticDet, Vdetected, Tdetectable, Tdetection, scanInterval, ";
    Ostream <<   "R, G, Mtotal, RecurrencesDetected, symptomaticDet, Vdetected, Tdetectable, Tdetection, scanInterval, ";
    Ostream <<   "R, G, Mtotal, RecurrencesDetected, symptomaticDet, Vdetected, Tdetectable, Tdetection, scanInterval, ";
    Ostream <<   "R, G, Mtotal, RecurrencesDetected, symptomaticDet, Vdetected, Tdetectable, Tdetection, scanInterval, " << endl;

    // calculate mean, average, min and max for all collumns
    string collumnNames[] = {"B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","AA","AB","AC","AD","AE","AF","AG","AH","AI","AJ","AK",};

    Ostream << "MEAN,";
    for (unsigned int c = 0; c < 36; c++) {
        Ostream << "=GEMIDDELDE(" << collumnNames[c] << "8:" << collumnNames[c] << 7+TOT << "),";
    }
    Ostream << endl << "STDEV,";
    for (unsigned int c = 0; c < 36; c++) {
        Ostream << "=STDEV(" << collumnNames[c] << "8:" << collumnNames[c] << 7+TOT << "),";
    }
    Ostream << endl << "MIN,";
    for (unsigned int c = 0; c < 36; c++) {
        Ostream << "=MIN(" << collumnNames[c] << "8:" << collumnNames[c] << 7+TOT << "),";
    }
    Ostream << endl << "MAX,";
    for (unsigned int c = 0; c < 36; c++) {
        Ostream << "=MAX(" << collumnNames[c] << "8:" << collumnNames[c] << 7+TOT << "),";
    }
    Ostream << endl << endl;

    return;
}

int GenerateNewPatient(int i, int j, string target)
{
    // *** MODEL PARAMETERS ***

    double Mu_Mtot = 18.5; // calibrated target 18.5 lower 23.5 upper 15.5
    double Sigma_Mtot = 15; // calibrated target 15 lower 16.5 upper 16.5
    int alfa_R = 9; // calibrated target 9 lower 12 upper 8
    int beta_R = 1; // calibrated (fixed on 1)
    double lambda_Symp = 0.00061; // calibrated 0.00061 lower 0.00042 upper 0.00081
    double lambda_Tdetectable = -0.00161; // calibrated target -0.00161 lower -0.00126 upper -0.00200
    double lambda_Growth = -0.0062; // derived from literature
    double VminDet = 0.06545; // derived from literature
    int oligoDefintion = 3; // literature
    double Sigma_FU = 15.5; // expert opinion (95% within 1 month of planned scan) (15.5)
    int minFUtime = 61; // expert opinion (2 months)
    int maxFUtime = 1825; // 5 years
    int FUschedule[] = {91, 182, 365, 548, 730, 1095, 1460, 1825}; // FU schedule based on guidelines

    double Sigma_Growth = 0; // sensitivity analysis parameter (10)
    double Simga_Detect = 0; // sensitivity analysis parameter (0.005)
    double Rho = 0; // correlation between Mtotal and G should be a value between 0 (no correlation) and 1 (maximum correlation) or -1 (negative?).

    // *** PATIENT SPECIFIC PARAMETERS ***

    int Mtotal = 0; // initialization
    while (Mtotal < 1) { Mtotal = ceil(rndm->Normal(Mu_Mtot,Sigma_Mtot)); } // estimate with 15.7% true oligo

    double R = 0; // initialization
    while (R == 0 || R == 1) { R = rndm->Quantile_Beta(alfa_R,beta_R); } // R cannot be 0 or 1. This should not happen, but does due to rounding errors.

    double G = lambda_Growth/log(((1-Rho)*rndm->Uniform(0.04,0.76)) + (Rho*(85.5-Mtotal)/112.5)); // *** adjusted for log(2) VDT G difference
    // both ranges(0.104-0.830) and limits VDT between 30 and 365 days
    double VDTrand = rndm->Normal(log(2)/G, Sigma_Growth); // random growth deviation of the smallest metastasis
    while (VDTrand < 30) { VDTrand = rndm->Normal(log(2)/G, Sigma_Growth); } // minimum of 30 days VDT

    int Tdetectable = 0;
    while (Tdetectable < 1 || Tdetectable > maxFUtime) { Tdetectable = ceil(log(rndm->Uniform(0.001,0.999))/lambda_Tdetectable); } // prevent log(0) from occurring

    int Tsymptomatic = maxFUtime; // initialization
    double randomMeta = log(rndm->Uniform(0.001, 0.999)); // prevent log(0) from occurring
    if (randomMeta > (1 - R*exp(lambda_Symp*Mtotal*Mtotal/(2*G)))) { Tsymptomatic = floor(Tdetectable + sqrt((2*randomMeta*log(R))/(lambda_Symp*G))); } // Mdetectable < Mtotal
    else { Tsymptomatic = floor(Tdetectable - Mtotal*log(R)/(2*G) - randomMeta/(lambda_Symp*Mtotal)); } // Mdetectable = Mtotal

    // *** loop through FU schedule to see if a scheduled scan happens before symptoms occur ***
    int Tdetection;
    int previousScan = 0;
    bool symptomaticDet = false;
    unsigned int fu = 0;
    int FUScanTime = floor(rndm->Normal(FUschedule[fu], Sigma_FU)); // FU scans are performed with some standard deviation from the guidelines
    while (FUScanTime < Tdetectable && fu != sizeof(FUschedule)) { // find first potential scan
        fu++;
        previousScan = FUScanTime; // store previous scan time here (it cannot be symptomatic)
        FUScanTime = floor(rndm->Normal(FUschedule[fu], Sigma_FU));
        if (FUScanTime < minFUtime) { FUScanTime = minFUtime; }
        if (FUScanTime > maxFUtime) { FUScanTime = maxFUtime; }
    }
    if(FUScanTime < Tsymptomatic) { Tdetection = FUScanTime; } else { symptomaticDet = true; Tdetection = Tsymptomatic; } // 14 days acceptable time not to plan unscheduled scan
    int RecurrencesDetected = 1 + floor((log(2)/VDTrand) *(Tdetectable-Tdetection)/log(R));
    if (RecurrencesDetected > Mtotal) { RecurrencesDetected = Mtotal; }
    double Vdetected = VminDet * pow(rndm->Normal(1,Simga_Detect),3) * exp((log(2)/VDTrand)*(Tdetection-Tdetectable)); // volume of largest detected tumor in cm3

    // *** write patient info to file (Simulation) ***
    if (target == "SIM") {
        // All Patients
        Ostream << ", " << R << ", " << G << ", " << Mtotal << ", " << RecurrencesDetected << ", " << symptomaticDet << ", " << Vdetected << ", " << Tdetectable << ", " << Tdetection << ", " << (Tdetection-previousScan) << ", ";

        // Poly Recurrence
        if (RecurrencesDetected > oligoDefintion) {
            Ostream << R << ", " << G << ", " << Mtotal << ", " << RecurrencesDetected << ", " << symptomaticDet << ", " << Vdetected << ", " << Tdetectable << ", " << Tdetection << ", " << (Tdetection-previousScan) << ", ";
        } else { Ostream << ",,,,,,,,,"; }

        // True Oligo Recurrence
        if (RecurrencesDetected <= oligoDefintion && Mtotal <= oligoDefintion) {
            Ostream << R << ", " << G << ", " << Mtotal << ", " << RecurrencesDetected << ", " << symptomaticDet << ", " << Vdetected << ", " << Tdetectable << ", " << Tdetection << ", " << (Tdetection-previousScan) << ", ";
        } else { Ostream << ",,,,,,,,,"; }

        // False Positive Oligo
        if (RecurrencesDetected <= oligoDefintion && Mtotal > oligoDefintion) {
            Ostream << R << ", " << G << ", " << Mtotal << ", " << RecurrencesDetected << ", " << symptomaticDet << ", " << Vdetected << ", " << Tdetectable << ", " << Tdetection << ", " << (Tdetection-previousScan) << ", ";
        } else { Ostream << ",,,,,,,,,"; }

        Ostream << endl;
    }

    // *** Calculate Calibration Target ***
    int calibrationTarget = 0;
    if (target == "MTOT" && Mtotal <= oligoDefintion) { calibrationTarget = 1; }
    if (target == "OLIGO" && RecurrencesDetected <= oligoDefintion) { calibrationTarget = 1; }
    if (target == "PFS") { calibrationTarget = floor(Tdetection/182.5); }
    if (target == "SYMP") {calibrationTarget = symptomaticDet; }

    return calibrationTarget;
}

void CloseOutputfile()
{
    Ostream.close();
    return;
}
