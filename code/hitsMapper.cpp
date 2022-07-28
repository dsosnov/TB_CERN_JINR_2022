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
    pair<string, string> run_pair = {"run_0091", "run47"};
    unsigned long long from = 0, to = 0;

    auto apvan = new apv(run_pair.second);
    apvan->useSyncSignal();
    auto hits_apv = apvan->GetCentralHits2ROnly(from, to);
    auto vmman = new evBuilder(run_pair.first, "g1_p25_s100-0&60", "map-20220721");
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
    out_APV.open("../out/APV_hits_all__.txt");
    out_VMM.open("../out/VMM_hits__.txt");
    int numOfMapped = 0;

    auto out = TFile::Open("../out/mapped_all.root", "recreate");

    auto stripsVMM = make_shared<TH1F>("stripsVMM", "stripsVMM", 360, 0, 360);
    auto mappedHitsPdo = make_shared<TH1F>("mappedHitsPdo", "mappedHitsPdo", 2000, 0, 2000);
    auto mappedHitsPdo_apv = make_shared<TH1F>("mappedHitsPdo_apv", "mappedHitsPdo_apv", 2000, 0, 2000);
    auto hitsPdo = make_shared<TH2F>("hitsPdo", "hitsPdo", 360, 0, 360, 2000, 0, 2000);
    auto hitsPdo_apv = make_shared<TH2F>("hitsPdo_apv", "hitsPdo_apv", 360, 0, 360, 2000, 0, 2000);
    auto hitsPdoAll = make_shared<TH1F>("hitsPdoAll", "hitsPdoAll", 2000, 0, 2000);
    auto hitsPdoAll_apv = make_shared<TH1F>("hitsPdoAll_apv", "hitsPdoAll_apv", 2000, 0, 2000);

    int a = -100;
    int b = -1;
    long long prevT = 0;
    for (unsigned long j = 0; j < hits_vmm_v.size(); j++)
    {
        if (hits_vmm_v.at(j).second.sync && a < 0)
        {   if (a == -100)
            {
                std::cout << "VMM first pule DAQ t = " << hits_vmm_v.at(j).second.timeFull() << "\n";
                prevT = hits_vmm_v.at(j).second.timeFull();
                a++;
                continue;
            }
            else 
            {
                std::cout << "VMM pulses delta T = " << hits_vmm_v.at(j).second.timeFull() - prevT << "\n";
                a++;
                prevT = hits_vmm_v.at(j).second.timeFull();
            }
        }
        if (hits_vmm_v.at(j).second.hitsX.size() != 0)
        {
            for (std::map<unsigned int, unsigned int>::iterator it = hits_vmm_v.at(j).second.hitsX.begin(); it != hits_vmm_v.at(j).second.hitsX.end(); ++it)
            {
                stripsVMM->Fill(it->first * (1 - 8e-3) - 8.46 / 0.25);
                int strip = it->first * (1 - 8e-3) - 8.46 / 0.25;
                int pdo = it->second;
                hitsPdo->Fill(strip, pdo);
                hitsPdoAll->Fill(pdo);
            }
        }
    }

    for (unsigned long i = 0; i < hits_apv_v.size(); i++)
    {
        if (hits_apv_v.at(i).second.sync && b < 0)
        {
            std::cout << "APV first pule DAQ t = " << hits_apv_v.at(i).second.timeFull() << "\n";
            b++;
        }
        if (hits_apv_v.at(i).second.hitsPerLayer.at(0).size() != 0)
        {
            for (auto &h : hits_apv_v.at(i).second.hitsPerLayer.at(0))
            {
                int strip = h.first * (1 - 8e-3) - 8.46 / 0.25;
                int pdo = h.second;
                hitsPdo_apv->Fill(strip, pdo);
                hitsPdoAll_apv->Fill(pdo);
            }
        }
    }
    // out->Write();
    // out->Close();

    // return 0;

    for (unsigned long i = 0; i < hits_apv_v.size() / 5; i++)
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
            std::cout << "Period " << nPeriodsAPV << "--- is sync! N = " << hits_apv_v.at(i).second.hitsPerLayer.at(0).size() << "\n";
        }
        else
        {
            if (prev_pulse_SRS != -1)
            {
                std::cout << "Not sync! N = " << hits_apv_v.at(i).second.hitsPerLayer.at(0).size() << "\n";
            }
        }

        if (prev_pulse_SRS != -1)
        {
            // std::cout << "\t Total:" << prev_pulse_SRS << "\t" << hits_apv_v.at(i).second.timeSrs << "\n";

            apv_hit_T = pulser_T + apv_time_from_SRS(prev_pulse_SRS, hits_apv_v.at(i).second.timeSrs);
            if (hits_apv_v.at(i).second.hitsPerLayer.at(0).size() != 0)
            {
                vector< pair <int,int> > apv_hits_vec;
                vector< pair <int,int> > apv_hits_vec_l0;
                vector< pair <int,int> > apv_hits_vec_l1;
                // std::cout << "---> with hits:" << prev_pulse_SRS << "\t" << hits_apv_v.at(i).second.timeSrs << "\n";
                // out_APV << "------- APV Period " << nPeriodsAPV << " -------- \n";
                for (auto &h: hits_apv_v.at(i).second.hitsPerLayer.at(2))
                {   
                    apv_hits_vec_l0.push_back(h);
                }
                for (auto &h: hits_apv_v.at(i).second.hitsPerLayer.at(1))
                {
                    int strip = h.first * (1 - 2.29e-3) - 2.412 / 0.25;
                    int pdo = h.second;
                    apv_hits_vec_l1.push_back(make_pair(strip, pdo));
                }
                for (auto &h: hits_apv_v.at(i).second.hitsPerLayer.at(0))
                {
                    int strip = h.first * (1 - 8e-3) - 8.46 / 0.25;
                    int pdo = h.second;
                    apv_hits_vec.push_back(make_pair(strip, pdo));
                }
                vector< pair <int,int> > vmm_hits_vec;

                int index = -200;
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
                            stripsVMM->Fill(it->first * (1 - 8e-3) - 8.46 / 0.25);
                            if (it->first == 171)
                                continue;
                            int strip = it->first * (1 - 8e-3) - 8.46 / 0.25;
                            int pdo = it->second;
                            vmm_hits_vec.push_back(make_pair(strip, pdo));
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
                            out_APV << "------- APV Double ReadOut Period " << nPeriodsAPV << " -------- T = " << T_apv - startT_pulse_apv << " [us] \n";
                            for (int k = 0; k < apv_hits_vec.size(); k++)
                            {
                                out_APV << k << "\t Strip: " << apv_hits_vec.at(k).first << "\n";
                                out_APV << "  \t PDO: " << apv_hits_vec.at(k).second << "\n";
                                mappedHitsPdo_apv->Fill(apv_hits_vec.at(k).second);
                            }
                            out_APV << "------- APV Layer 0 -------- \n";
                            for (int k = 0; k < apv_hits_vec_l0.size(); k++)
                            {
                                out_APV << k << "\t Strip: " << apv_hits_vec_l0.at(k).first << "\n";
                                out_APV << "  \t PDO: " << apv_hits_vec_l0.at(k).second << "\n";
                            }
                            out_APV << "------- APV Layer 1 -------- \n";
                            for (int k = 0; k < apv_hits_vec_l1.size(); k++)
                            {
                                out_APV << k << "\t Strip: " << apv_hits_vec_l1.at(k).first << "\n";
                                out_APV << "  \t PDO: " << apv_hits_vec_l1.at(k).second << "\n";
                            }

                            out_APV << "------- VMM Period " << nPeriods / 200 << "  (" << nPeriods % 200 << ") -------- T = " << pulseTime + diff_hit * 25 / 1000 << "[us] (dT = " << dt_apv_vmm << " [us]) \n";

                            for (int l = 0; l < vmm_hits_vec.size(); l++)
                            {
                                out_APV << "\t Strip: " << vmm_hits_vec.at(l).first << "\n";
                                out_APV << "\t PDO: " << vmm_hits_vec.at(l).second << "\n";
                                mappedHitsPdo->Fill(vmm_hits_vec.at(l).second);
                            }

                        }
                    }
                }
                
            }
        }
    }
    out_APV.close();
    out_VMM.close();
    out->Write();
    out->Close();
}
