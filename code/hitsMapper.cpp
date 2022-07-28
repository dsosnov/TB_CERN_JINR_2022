#include "apv.C"
#include "evBuilder.C"
#include <iostream>
#include <fstream>


int apv_time_from_SRS(int srs1, int srs2)
{
    int diff = 0;
    int srsT_period = 16777215;
    if (srs2 - srs1 >= 0)
    {
        diff = srs2 - srs1;
    }
    else
    {
        diff = srs2 + srsT_period - srs1;
    }

    return diff * 25 / 1000;
}

void hitsMapper()
{
    pair<string, string> run_pair = {"run_0832_cut", "run423_cut"};
    unsigned long long from = 0, to = 0;

    auto apvan = new apv(run_pair.second);
    apvan->useSyncSignal();
    auto hits_apv = apvan->GetCentralHits2ROnly(from, to);
    auto vmman = new evBuilder(run_pair.first, "g1_p25_s100-0&60", "map-20220605");
    vmman->useSyncSignal();
    auto hits_vmm = vmman->GetCentralHits(from, to);

    vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> hits_vmm_v;
    hits_vmm_v.assign(hits_vmm.begin(), hits_vmm.end());
    vector<pair<unsigned long, apv::doubleReadoutHits>> hits_apv_v;
    hits_apv_v.assign(hits_apv.begin(), hits_apv.end());

    long long startT_apv = 0;
    long long startT_pulse_apv = 0;
    long long T_apv = 0;
    int nPeriodsAPV = 0;

    int prevSRS = -1;
    int prev_pulse_SRS = -1;
    long long pulser_T = 0;
    long long apv_hit_T = 0;

    ofstream out_APV;
    ofstream out_VMM;
    out_APV.open("../out/APV_hits__.txt");
    out_VMM.open("../out/VMM_hits__.txt");
    int numOfMapped = 0;

    


    for (unsigned long i = 0; i < hits_apv_v.size(); i++)
    {
        if (prevSRS == -1)
        {
            prevSRS = hits_apv_v.at(i).second.timeSrs;
            startT_apv = 25 * hits_apv_v.at(i).second.timeSrs / 1000;
            T_apv = startT_apv;
        }
        else
        {
            T_apv += apv_time_from_SRS(prevSRS, hits_apv_v.at(i).second.timeSrs);
            prevSRS = hits_apv_v.at(i).second.timeSrs;
        }

        if (hits_apv_v.at(i).second.sync)
        {
            if (prev_pulse_SRS == -1)
            {
                startT_pulse_apv = T_apv;
                prev_pulse_SRS = hits_apv_v.at(i).second.timeSrs;
            }
            else
            {
                nPeriodsAPV = round((T_apv - startT_pulse_apv) * 1.0 / 1e4);
                prev_pulse_SRS = hits_apv_v.at(i).second.timeSrs;
                pulser_T = nPeriodsAPV * 1e4;
            }
            std::cout << "Period " << nPeriodsAPV << "--- is sync! N = " << hits_apv_v.at(i).second.hits.size() << "\n";
        }
        else
        {
            if (prev_pulse_SRS != -1)
            {
                std::cout << "Not sync! N = " << hits_apv_v.at(i).second.hits.size() << "\n";
            }
        }

        if (prev_pulse_SRS != -1)
        {
            // std::cout << "\t Total:" << prev_pulse_SRS << "\t" << hits_apv_v.at(i).second.timeSrs << "\n";

            apv_hit_T = pulser_T + apv_time_from_SRS(prev_pulse_SRS, hits_apv_v.at(i).second.timeSrs);
            if (hits_apv_v.at(i).second.hits.size() != 0)
            {
                vector< pair <int,int> > apv_hits_vec;
                // std::cout << "---> with hits:" << prev_pulse_SRS << "\t" << hits_apv_v.at(i).second.timeSrs << "\n";
                // out_APV << "------- APV Period " << nPeriodsAPV << " -------- \n";
                for (int j = 0; j < hits_apv_v.at(i).second.hits.size(); j++)
                {
                    apv_hits_vec.push_back(make_pair(hits_apv_v.at(i).second.hits.at(j).strip, hits_apv_v.at(i).second.hits.at(j).max_q));
                    // out_APV << j << "\t Strip: " << hits_apv_v.at(i).second.hits.at(j).strip << "\n";
                    // out_APV << "  \t PDO: " << hits_apv_v.at(i).second.hits.at(j).max_q << "\n";
                    // out_APV << "  \t hit T: " << T_apv - startT_pulse_apv << "\n";
                }

                vector< pair <int,int> > vmm_hits_vec;

                int index = -652;
                int bad = 0;
                int prevSyncBcid = 0;
                int prevPrevSyncBcid = 0;
                int nPeriods = 0;

                long long pulseTime = 0;
                for (unsigned long j = 0; j < hits_vmm_v.size(); j++)
                {

                    if (hits_vmm_v.at(j).second.sync)
                    {
                        if (index < 0)
                        {
                            index++;
                            continue;
                        }
                        if (index == 0)
                        {
                            prevSyncBcid = hits_vmm_v.at(j).second.bcid;
                            index++;
                        }
                        else
                        {
                            int diff = 0;
                            if (hits_vmm_v.at(j).second.bcid - prevSyncBcid > 0)
                            {
                                diff = hits_vmm_v.at(j).second.bcid - prevSyncBcid;
                            }
                            else
                            {
                                diff = hits_vmm_v.at(j).second.bcid + 4096 - prevSyncBcid;
                            }
                            int diffDiff = 0;

                            if (hits_vmm_v.at(j).second.bcid - prevPrevSyncBcid > 0)
                            {
                                diffDiff = hits_vmm_v.at(j).second.bcid - prevPrevSyncBcid;
                            }
                            else
                            {
                                diffDiff = hits_vmm_v.at(j).second.bcid + 4096 - prevPrevSyncBcid;
                            }

                            if (!(diff >= 2000 - 5 && diff <= 2000 + 5) && !(diff >= 4000 - 5 && diff <= 4000 + 5))
                            {
                                // if (diff == 1904)
                                // {
                                //     diff = 6000;
                                // }
                                // else if (diff == 3904)
                                // {
                                //     diff = 8000;
                                // }
                                // else 
                                if (!(diffDiff >= 2000 - 5 && diffDiff <= 2000 + 5) && !(diffDiff >= 4000 - 5 && diffDiff <= 4000 + 5))
                                {
                                    bad++;
                                    prevPrevSyncBcid = prevSyncBcid;
                                    prevSyncBcid = hits_vmm_v.at(j).second.bcid;
                                    continue;
                                }
                                else
                                {
                                    diff = diffDiff;
                                }
                            }

                            prevPrevSyncBcid = prevSyncBcid;
                            prevSyncBcid = hits_vmm_v.at(j).second.bcid;
                            nPeriods += round(diff * 1.0 / 2000.0);
                            pulseTime = nPeriods * 50;
                        }

                    }
                    

                    if (nPeriods / 200 == nPeriodsAPV && index > 0 && hits_vmm_v.at(j).second.hitsX.size() != 0)
                    {
                        int diff_hit = 0;
                        if (hits_vmm_v.at(j).second.bcid - prevSyncBcid >= 0)
                        {
                            diff_hit = hits_vmm_v.at(j).second.bcid - prevSyncBcid;
                        }
                        else
                        {
                            diff_hit = hits_vmm_v.at(j).second.bcid + 4096 - prevSyncBcid;
                        }
                        int dt_apv_vmm = T_apv - startT_pulse_apv - pulseTime;
                        // std::cout << nPeriods / 200 << " \t " << hits_vmm_v.at(j).second.hitsX.size() << "\n";
                        // out_APV << "------- VMM Period " << nPeriods / 200 << "  (" << nPeriods % 200 << ") -------- dT = " << dt_apv_vmm << "\n";
                        if (abs(dt_apv_vmm) > 800)
                            continue;

                        dt_apv_vmm -= diff_hit * 25 / 1000;
                        for(std::map<unsigned int,unsigned int>::iterator it = hits_vmm_v.at(j).second.hitsX.begin(); it != hits_vmm_v.at(j).second.hitsX.end(); ++it){
                            if (it->first == 178 || it->first == 186)
                                continue;
                            vmm_hits_vec.push_back(make_pair(it->first, it->second));
                            // out_APV << "\t Strip: " << it->first << "\n";
                            // out_APV << "\t PDO: " << it->second << "\n";
                        }
                        bool hitMapped = false;

                        for (int k = 0; k < apv_hits_vec.size(); k++)
                        {
                            if (hitMapped)
                                continue;

                            if(vmm_hits_vec.size() != 0 && abs(apv_hits_vec.at(k).first - vmm_hits_vec.at(0).first) < 5 && vmm_hits_vec.at(0).first)
                                hitMapped = true;
                        }

                        if (hitMapped)
                        {
                            numOfMapped++;
                            std::cout << "!!! Mapped hits: " << numOfMapped << " !!! \n";
                            out_APV << "------- APV Period " << nPeriodsAPV << " -------- \n";
                            for (int k = 0; k < apv_hits_vec.size(); k++)
                            {
                                out_APV << k << "\t Strip: " << apv_hits_vec.at(k).first << "\n";
                                out_APV << "  \t PDO: " << apv_hits_vec.at(k).second << "\n";
                            }

                            out_APV << "------- VMM Period " << nPeriods / 200 << "  (" << nPeriods % 200 << ") -------- dT = " << dt_apv_vmm << "\n";

                            for (int l = 0; l < vmm_hits_vec.size(); l++)
                            {
                                out_APV << "\t Strip: " << vmm_hits_vec.at(l).first << "\n";
                                out_APV << "\t PDO: " << vmm_hits_vec.at(l).second << "\n";
                            }

                        }
                    }
                }
                
            }
        }
    }
    out_APV.close();
    out_VMM.close();
}
