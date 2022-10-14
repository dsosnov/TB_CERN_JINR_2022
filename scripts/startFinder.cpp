#include "../code/apv.C"
#include "../code/evBuilder.C"
#include <iostream>
#include <fstream>

#include <vector>
#include <utility>
#include <tuple>
#include <map>
#include <set>
#include <optional>

using std::cout, std::endl;
using std::ifstream, std::ofstream;

using std::vector;
using std::pair, std::make_pair;
using std::tuple, std::get;
using std::map;
using std::set;
using std::optional, std::nullopt;

optional<int> calculateVMMNPulsers(int bcidDiff, int maxDiff = 2, int maxNPulsers = 9)
{
    bcidDiff += (bcidDiff > 0) ? 0 : 4096;
    int predicted;
    for(auto i = 1; i <= maxNPulsers; i++)
    {
      predicted = (2000 * i) % 4096;
      if(abs(predicted - bcidDiff) <= maxDiff)
        return i;
    }
    return nullopt;
}

void startFinder(string runVMM = "0832", string runAPV = "423", int firstApvPulse = 0)
{
    pair<string, string> run_pair = {Form("run_%s", runVMM.c_str()), Form("run%s", runAPV.c_str())};

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
    int pulseN = -firstApvPulse;

    for (unsigned long i = 0; i < apvan->GetEntries(); i++)
    {
        auto hit_apv = apvan->GetCentralHits2ROnlyData(i);

        if (hit_apv.sync)
        {
            T_apv = hit_apv.timeSec*1E6 + hit_apv.timeMSec;
            if (pulseN == 0)
            {
                first_pulse_T_apv = T_apv;
            }
            cout << "APV Pulse event " << i << " N = " << pulseN << "\t T = " << T_apv << "\t from last pulse dT = " << ((prev_T_apv > 0) ? T_apv - prev_T_apv : 0) << endl;
            prev_T_apv = T_apv;
            if (pulseN == 10)
                break;
            pulseN++;
        }
    }

    cout << "\n\n";

    pulseN = 0;

    constexpr bool checkInTimeRange = true; // true - first in range 300 ms from APV; false - select first events;
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

            if(checkInTimeRange)
            {
                if (abs(T_vmm - first_pulse_T_apv) < 1000)
                {
                    cout << "VMM Pulse event " << i << " N = " << pulseN << "\t T = " << T_vmm << "\t from last pulse dT = " << T_vmm - prev_T_vmm<< "\t vs FirstAPV = " << T_vmm - first_pulse_T_apv << endl;
                } else if(T_vmm > first_pulse_T_apv)
                {
                    break;
                }
            }
            else
            {
                cout << "VMM Pulse event " << i << " N = " << pulseN << "\t T = " << T_vmm << "\t from last pulse dT = " << T_vmm - prev_T_vmm<< "\t vs FirstAPV = " << T_vmm - first_pulse_T_apv << endl;
                if (pulseN >= 20 || T_vmm > first_pulse_T_apv + 300) // 100000
                  break;
            }
            prev_T_vmm = T_vmm;
            pulseN++;
        }
    }

    cout << "First APV pulse T = " << first_pulse_T_apv << endl;
    cout << "First VMM pulse T = " << first_pulse_T_vmm << endl;

    cout << "First APV vs VMM dt = " << first_pulse_T_apv - first_pulse_T_vmm << " or " << (first_pulse_T_apv - first_pulse_T_vmm) / 10000 << " pulses" << endl;


}
