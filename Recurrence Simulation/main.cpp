#include "simulation.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    cout << endl
    << "This Non-Small Cell Lung Cancer micro-simulation model was made by" << endl
    << "by Harold Wolff and Veerle Coupe" << endl
    << "VU University Medical Center" << endl
    << "The Netherlands (2018)" << endl << endl;

    // *** loop through parameters if desired, otherwise set to i < 1 && j < 1 ***

    // either run a simulation or calibrate the model
    string target = "SIM"; // SIM SYMP MTOT OLIGO PFS
    double MSE = 0;
    unsigned int maxPatients = 100000;
    unsigned int I_max = 2; // loop size
    unsigned int J_max = 2; // loop size also used for calibration

    if (target == "SIM") { cout << "generating new output files ..." << endl; }
    else { cout << "calculating MSE..." << endl; }

    for (unsigned int i = 1; i < I_max; i++) {
        for (unsigned int j = 1; j < J_max; j++) {

            // add rules to calculate I and J from i and j here.
            int I = i;
            int J = j;

            if (target == "SIM") { GenerateOutputfile(maxPatients, I, J); }

            // reset calibration targets
            double simVAR = 0;
            double SSE1 = 0;
            double SSE2 = 0;
            double SSE3 = 0;
            double SSE4 = 0;
            double SSE5 = 0;
            double SSE6 = 0;
            double SSE7 = 0;
            double SSE8 = 0;
            double SSE9 = 0;

            // run simulation
            for (unsigned int patient = 0; patient < maxPatients; patient++) {
                if (target == "PFS") {
                    simVAR = GenerateNewPatient(I, J, target);
                    if (simVAR > 0) { SSE1 +=1; }
                    if (simVAR > 1) { SSE2 +=1; }
                    if (simVAR > 2) { SSE3 +=1; }
                    if (simVAR > 3) { SSE4 +=1; }
                    if (simVAR > 4) { SSE5 +=1; }
                    if (simVAR > 5) { SSE6 +=1; }
                    if (simVAR > 6) { SSE7 +=1; }
                    if (simVAR > 7) { SSE8 +=1; }
                    if (simVAR > 8) { SSE9 +=1; }
                } else { simVAR += GenerateNewPatient(I, J, target); }
            }

            if (target == "SIM") { CloseOutputfile(); }
            if (target == "MTOT") { MSE += pow(((simVAR/maxPatients) - 0.0471), 2); } // oligos without micrometstases target 0.0471 lower 0.0330 upper 0.0611
            if (target == "OLIGO") { MSE += pow(((simVAR/maxPatients) - 0.2999), 2); } // oligos with and without micrometastases target 0.2999 lower 0.2758 upper 0.3239
            if (target == "SYMP") { MSE += pow(((simVAR/maxPatients) - 0.3350), 2);} // symptomatic detected recurrences target 0.3350 lower 0.2884 upper 0.3817
            if (target == "PFS") {
                // PFS at specific time points and weights for patients at risk
                MSE += (pow(((SSE1/maxPatients) - 0.82751), 2) * 0.214); // target 0.82751 lower 0.79586 upper 0.85917
                MSE += (pow(((SSE2/maxPatients) - 0.56114), 2) * 0.183); // target 0.56114 lower 0.51261 upper 0.60967
                MSE += (pow(((SSE3/maxPatients) - 0.43502), 2) * 0.162); // target 0.43502 lower 0.38110 upper 0.48893
                MSE += (pow(((SSE4/maxPatients) - 0.30767), 2) * 0.132); // target 0.30767 lower 0.24866 upper 0.36667
                MSE += (pow(((SSE5/maxPatients) - 0.25452), 2) * 0.105); // target 0.25452 lower 0.19294 upper 0.31610
                MSE += (pow(((SSE6/maxPatients) - 0.19142), 2) * 0.083); // target 0.19142 lower 0.12629 upper 0.25655
                MSE += (pow(((SSE7/maxPatients) - 0.16143), 2) * 0.061); // target 0.16143 lower 0.09372 upper 0.22913
                MSE += (pow(((SSE8/maxPatients) - 0.00000), 2) * 0.031); // target 0.00000 lower 0.00000 upper 0.08991
                MSE += (pow(((SSE9/maxPatients) - 0.00000), 2) * 0.031); // target 0.00000 lower 0.00000 upper 0.08991
            }
        }
    }
    if (target == "SIM") { cout << "simulations finished ..." << endl; }
    else { cout << "MSE: " << MSE << endl; }

    return 0;
}
