#include "apv.C"
#include "evBuilder.C"
#include <iostream>
#include <fstream>

#include <vector>
#include <utility>
#include <tuple>
#include <map>
#include <set>

using std::cout, std::endl;
using std::ifstream, std::ofstream;

using std::vector;
using std::pair, std::make_pair;
using std::tuple, std::get;
using std::map;
using std::set;


void startFinder()
{
    pair<string, string> run_pair = {"run_0087", "run43"};

    auto apvan = new apv(run_pair.second);
    apvan->useSyncSignal();
    auto vmman = new evBuilder(run_pair.first, "g1_p25_s100-0&60", "map-20220605");
    vmman->useSyncSignal();

    cout << "Num of events: APV -- " << apvan->GetEntries() << "; VMM -- " << vmman->GetEntries() << endl;

    unsigned long i;
    long long first_pulse_T_apv = 0;
    long long first_pulse_T_vmm = 0;
    long long T_apv = 0;
    long long T_vmm = 0;
    long long prev_T_apv = 0;
    long long prev_T_vmm = 0;
    int pulseN = 0;

    for (unsigned long i = 0; i < apvan->GetEntries(); i++)
    {
        auto hit_apv = apvan->GetCentralHits2ROnlyData(i);

        if (hit_apv.sync)
        {
            if (first_pulse_T_apv == 0)
            {
                first_pulse_T_apv = hit_apv.timeSec*1E6 + hit_apv.timeMSec;
                prev_T_apv = first_pulse_T_apv;
            }
            T_apv = hit_apv.timeSec*1E6 + hit_apv.timeMSec;

            cout << "APV Pulse event " << i << " N = " << pulseN << "\t T = " << T_apv << "\t from last pulse dT = " << T_apv - prev_T_apv << endl;
            prev_T_apv = T_apv;
            if (pulseN == 20)
                break;
            pulseN++;
        }
    }

    cout << "\n\n";

    pulseN = 0;

    for (unsigned long i = 0; i < vmman->GetEntries(); i++)
    {
        auto hit_vmm = vmman->GetCentralHitsData(i);

        if (hit_vmm.sync)
        {
            if (first_pulse_T_vmm == 0)
            {
                first_pulse_T_vmm = hit_vmm.timeSec*1E6 + hit_vmm.timeMSec;
                prev_T_vmm = first_pulse_T_vmm;
            }
            T_vmm = hit_vmm.timeSec*1E6 + hit_vmm.timeMSec;

            // if (abs(T_vmm - first_pulse_T_apv) < 300)
            // {
            //     cout << "VMM Pulse event " << i << " N = " << pulseN << "\t T = " << T_vmm << "\t from last pulse dT = " << T_vmm - prev_T_vmm<< "\t vs FirstAPV = " << T_vmm - first_pulse_T_apv << endl;
            // }
            cout << "VMM Pulse event " << i << " N = " << pulseN << "\t T = " << T_vmm << "\t from last pulse dT = " << T_vmm - prev_T_vmm<< "\t vs FirstAPV = " << T_vmm - first_pulse_T_apv << endl;
            prev_T_vmm = T_vmm;
            if (pulseN == 10000)
                break;
            pulseN++;
        }
    }

    cout << "First APV pulse T = " << first_pulse_T_apv << endl;
    cout << "First VMM pulse T = " << first_pulse_T_vmm << endl;

    cout << "First APV vs VMM dt = " << first_pulse_T_apv - first_pulse_T_vmm << " or " << (first_pulse_T_apv - first_pulse_T_vmm) / 10000 << " pulses" << endl;


}