#include<omp.h>
#include "ehhpop.h"
#include "CmdLine.h"

using namespace std;

int main(int argc, char *argv[])
{
    bool distcorrect = true;
    CCmdLine cmdline;
    string pop1in, pop2in, mapin;
    if (cmdline.SplitLine(argc, argv) < 1) {
        // no switches were given on the command line, abort
        cout << "No input files\n";
        cout << "-h [input file 1] [input file 2]\n";
        cout << "-m [map file]\n";
        exit(-1);
    }
    if (!cmdline.HasSwitch("-h")) {
        cout << "ERROR: no haplotype files\n";
        cout << "-h [input file 1] [input file 2]\n";
        cout << "-m [map file]\n";
        exit(1);
    } else {
        pop1in = cmdline.GetArgument("-h", 0);
        pop2in = cmdline.GetArgument("-h", 1);
    }
    if (!cmdline.HasSwitch("-m")) {
        cout << "ERROR: no map file\n";
        cout << "-h [input file 1] [input file 2]\n";
        cout << "-m [map file]\n";
        exit(1);
    } else{
        mapin = cmdline.GetArgument("-m", 0);
    }

    if (cmdline.HasSwitch("-nd"))
        distcorrect = false;

    Ehhpop pop1(pop1in.c_str(), mapin.c_str());
    Ehhpop pop2(pop2in.c_str(), mapin.c_str());
    Ehhpop combined(pop1, pop2);

    // OK, the OMP part, we want to declare all variables outide of the loop.
    string rs;
    vector<string> reps(combined.Nsnp);
    char s[512], stmp[512];
    int i, pos, left, right;
    double p1_left, p1_right, IA, p2_left, p2_right, IB, logratio;

    #pragma omp parallel shared(distcorrect, combined, pop1, pop2, reps) private(i, rs, pos, left, right, p1_left, p1_right, IA, p2_left, p2_right, IB, logratio, s, stmp)
    {
        #pragma omp for
        for(i = 0; i<combined.Nsnp; i++) {
            rs = combined.index2rs[i];
            pos = combined.rs2pos[rs];
            left = combined.findcutoff(i, 0.05, false, distcorrect);
            right = combined.findcutoff(i, 0.05, true, distcorrect);
            //cout << left << " "<< right << "\n";
            if(left==-1 || right ==-1) {
                continue;
            }

            p1_left = pop1.integrate_ehh(i, left, false, distcorrect);
            p1_right = pop1.integrate_ehh(i, right, true, distcorrect);
            // cout << p1_left << " "<< p1_right << " p1\n";
            sprintf(s, "%4.6f %4.6f p1\n", p1_left, p1_right);
            IA = p1_left+p1_right;

            p2_left = pop2.integrate_ehh(i, left, false, distcorrect);
            p2_right = pop2.integrate_ehh(i, right, true, distcorrect);
            IB = p2_left+p2_right;
            // cout << p2_left << " "<< p2_right << " p2\n";
            sprintf(stmp, "%4.6f %4.6f p2\n", p2_left, p2_right);
            strcat(s, stmp);
            logratio = log(IA/IB);
            // cout << rs << " " << pos << " "<< IA << " " << IB << " " << logratio << "\n";
            sprintf(stmp, "%s %i %4.6f %4.6f %4.6f\n", rs.c_str(), pos, IA, IB, logratio);
            strcat(s, stmp);
            reps[i] = s;
        }
    } // omp parallel region closes.

    for(i = 0; i<combined.Nsnp; i++)
        cout << reps[i];

    return 0;
}
