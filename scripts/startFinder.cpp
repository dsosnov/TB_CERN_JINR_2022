#include "../code/apv.C"
#include "../code/evBuilder.C"
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
    return -1;
}



void startFinder()
{
    pair<string, string> run_pair = {"run_0080", "run31"};

    auto apvan = new apv(run_pair.second);
    apvan->useSyncSignal();
    auto vmman = new evBuilder(run_pair.first, "g1_p25_s100-0&60", "map-20220721");
    vmman->useSyncSignal();

    cout << "Num of events: APV -- " << apvan->GetEntries() << "; VMM -- " << vmman->GetEntries() << endl;

    unsigned long i;
    long long first_pulse_T_apv = 0;
    long long first_pulse_T_vmm = 0;
    long long T_apv = 0;
    long long T_vmm = 0;    
    long long T_vmm_pred = 0;
    int bcid_vmm = 0;
    long long prev_T_apv = 0;
    long long prev_T_vmm = 0;
    int prev_bcid_vmm = 0;
    int pulseN = 0;
    int pulseNVMM = 0;
    int srsTs = 0;
    int prevSrsrTs = 0;


    for (unsigned long i = 0; i < apvan->GetEntries(); i++)
    {
        auto hit_apv = apvan->GetCentralHits2ROnlyData(i);

        if (hit_apv.sync)
        {
            srsTs = hit_apv.timeSrs;
            if (first_pulse_T_apv == 0)
            {
                first_pulse_T_apv = hit_apv.timeSec*1E6 + hit_apv.timeMSec;
                prev_T_apv = first_pulse_T_apv;
                prevSrsrTs = hit_apv.timeSrs;
            }
            T_apv = hit_apv.timeSec*1E6 + hit_apv.timeMSec;

            if (true)
            {
                cout << "APV Pulse event " << i << " N = " << pulseN << "\t T = " << (T_apv - first_pulse_T_apv)/1e6 << "\t from last pulse dT = " << T_apv - prev_T_apv << endl;
                cout << "delta SRS TS = " << srsTs - prevSrsrTs << "\n";
            }

            prev_T_apv = T_apv;
            prevSrsrTs = srsTs;
            if (pulseN > 20)
                break;
            pulseN++;
        }
    }

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
            }
            T_vmm = hit_vmm.timeSec*1E6 + hit_vmm.timeMSec;
            int diff = (bcid_vmm - prev_bcid_vmm > 0) ? bcid_vmm - prev_bcid_vmm : bcid_vmm + 4096 - prev_bcid_vmm;
            auto npulsers = calculateVMMNPulsers(diff, 1, 86);
            T_vmm_pred +=  npulsers * 50;

            // 100.827 dt_apv_vmm = 1658668867312046 - 1658668857229304 - 10082950 - round(0 * 25.0 / 1000.0) = -379
            // 100.884 dt_apv_vmm = 1658668867317731 - 1658668857229304 - 10088900 - round(0 * 25.0 / 1000.0) = -644
            // 109.177 dt_apv_vmm = 1658668868146986 - 1658668857229304 - 10918500 - round(0 * 25.0 / 1000.0) = -988

            // VMM Pulse event 395025 N = 161217        T = 8.57976     from last pulse dT (DAQ) = 186613 | dT (BCID) = 50      N of hits inside 5873
            // DAQ vs BCID dt = 187010  delta BCID = 2000

            // if (npulsers < 0 || hit_vmm.pdo != 1012)
            //     continue;

            // if ((T_vmm - first_pulse_T_vmm)/1e6 > 0)
            if (true) 
            {
                // cout << "VMM Pulse event " << i << " N = " << pulseN << "\t T = " << T_vmm << "\t from last pulse dT = " << T_vmm - prev_T_vmm<< "\t vs FirstAPV = " << T_vmm - first_pulse_T_apv << endl;
                cout << "VMM Pulse event " << i << " N = " << pulseN << "\t T = " << (T_vmm - first_pulse_T_vmm)/1e6 << "\t from last pulse dT (DAQ) = " << T_vmm - prev_T_vmm << " | dT (BCID) = " << npulsers * 50 << "\t N of hits inside " << numOfVmmHits << endl;
                cout << "DAQ vs BCID dt = " << T_vmm - first_pulse_T_vmm - T_vmm_pred << "\t delta BCID = " << diff <<  "\t vs FirstAPV = " << T_vmm - first_pulse_T_apv << "\n";

            }
            // cout << "VMM Pulse event " << i << " N = " << pulseN << "\t T = " << T_vmm << "\t from last pulse dT = " << T_vmm - prev_T_vmm<< "\t vs FirstAPV = " << T_vmm - first_pulse_T_apv << endl;
            prev_T_vmm = T_vmm;
            prev_bcid_vmm = bcid_vmm;
            if (pulseN > 300)
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