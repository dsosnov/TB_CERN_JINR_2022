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

int calculateVMMNPulsers(int bcidDiff, int maxDiff = 2, int maxNPulsers = 9)
{
    bcidDiff += (bcidDiff > 0) ? 0 : 4096;
    int predicted;
    for(auto i = 1; i <= maxNPulsers; i++)
    {
      predicted = (2000 * i) % 4096;
      if(abs(predicted - bcidDiff) <= maxDiff)
        return i;
    }
    return 0;
}

void rawAna()
{
    pair<string, string> run_pair = {"run_0826", "run418"};

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
    long long T_vmm_pred = 0;
    long long bcid_vmm = 0;
    long long prev_T_apv = 0;
    long long prev_T_vmm = 0;
    long long prev_bcid_vmm = 0;
    int pulseN = 0;
    int pulseNVMM = 0;
    bool checkInTimeRange = true; // true - first in range 300 ms from APV; false - select first events;

    // for (unsigned long i = 0; i < apvan->GetEntries(); i++)
    // {
    //     auto hit_apv = apvan->GetCentralHits2ROnlyData(i);

    //     if (hit_apv.sync)
    //     {
    //         if (first_pulse_T_apv == 0)
    //         {
    //             first_pulse_T_apv = hit_apv.timeSec*1E6 + hit_apv.timeMSec;
    //             prev_T_apv = first_pulse_T_apv;
    //         }
    //         T_apv = hit_apv.timeSec*1E6 + hit_apv.timeMSec;
    //         if (T_apv - prev_T_apv > 0 && T_apv - first_pulse_T_apv > 359e6)
    //         {
    //             cout << "APV Pulse event " << i << " N = " << pulseN << "\t T = " << (T_apv - first_pulse_T_apv)/1e6 << "\t from last pulse dT = " << (T_apv - prev_T_apv)/1e6  << endl;
    //         }
    //         prev_T_apv = T_apv;
    //         // if (pulseN == 10)
    //         if (T_apv - first_pulse_T_apv > 380e6)
    //             break;
    //         pulseN++;
    //     }
    // }

    cout << "\n\n";

    pulseN = 0;
    
    int numOfVmmHits = 0;

    for (unsigned long i = 0; i < vmman->GetEntries(); i++)
    {
        auto hit_vmm = vmman->GetCentralHitsData(i);

        if (hit_vmm.sync)
        {
            bcid_vmm = hit_vmm.bcid;

            numOfVmmHits += hit_vmm.hitsX.size();
            
            if (first_pulse_T_vmm == 0)
            {
                first_pulse_T_vmm = hit_vmm.timeSec*1E6 + hit_vmm.timeMSec;
                prev_T_vmm = first_pulse_T_vmm;
                prev_bcid_vmm = bcid_vmm;
            }

            T_vmm = hit_vmm.timeSec*1E6 + hit_vmm.timeMSec;
            int diff = (bcid_vmm - prev_bcid_vmm > 0) ? bcid_vmm - prev_bcid_vmm : bcid_vmm + 4096 - prev_bcid_vmm;
            auto npulsers = calculateVMMNPulsers(diff, 1, 30);
            if (npulsers < 0 || hit_vmm.pdo > 966 || hit_vmm.pdo < 947)
            {
                // hbcidDiffIgnored->Fill(diff);
                continue;
            }
            T_vmm_pred += npulsers * 50;


            if (T_vmm - prev_T_vmm > 200 && T_vmm - first_pulse_T_vmm > 359e6)
            {
                cout << "VMM Pulse event " << i << " N = " << pulseN << "\t T = " << (T_vmm - first_pulse_T_vmm)/1e6 << "\t from last pulse dT (DAQ) = " << T_vmm - prev_T_vmm << " | dT (BCID) = " << npulsers * 50 << "\t N of hits inside " << numOfVmmHits << endl;
                cout << "DAQ vs BCID dt = " << T_vmm - first_pulse_T_vmm - T_vmm_pred << "\t delta BCID = " << diff << "\n";
            }
            // numOfVmmHits = 0;
            prev_T_vmm = T_vmm;
            prev_bcid_vmm = bcid_vmm;
            // if (pulseN == 1000)
            if (T_vmm - first_pulse_T_vmm > 380e6)
                break;
            pulseN++;
        }
        else
        {
            numOfVmmHits += hit_vmm.hitsX.size();
        }
    }

    cout << "First APV pulse T = " << first_pulse_T_apv << endl;
    cout << "First VMM pulse T = " << first_pulse_T_vmm << endl;

    cout << "First APV vs VMM dt = " << first_pulse_T_apv - first_pulse_T_vmm << " or " << (first_pulse_T_apv - first_pulse_T_vmm) / 10000 << " pulses" << endl;


}