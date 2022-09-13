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

    return round(diff * 25.0 / 1000.0);
}

void hitsMapper_reverse()
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

    ofstream out_txt;
    out_txt.open("../out/hits_maped_0832_423.txt");

    auto out = TFile::Open("../out/reverse_mapped_0832_423.root", "recreate");

    auto stripsVMM = make_shared<TH1F>("stripsVMM", "stripsVMM", 360, 0, 360);
    auto mappedHitsPdo = make_shared<TH1F>("mappedHitsPdo", "mappedHitsPdo", 2000, 0, 2000);

    int index = 0;
    int prevSyncBcid = 0;
    int prevPrevSyncBcid = 0;
    int nPeriods = 0;

    long long pulseTime = 0;
    long long vmmHitTime = 0;
    vector<pair<int, int>> vmm_hits_vec;
    int numOfMapped = 0;

    double vmmPrevHitT_mapped = 0;
    double apvPrevHitT_mapped = 0;

    int n_of_unmapped_becauese_of_T = 0;
    int num = 0;
    for (unsigned long j = 0; j < hits_vmm_v.size(); j++)
    {
        if (hits_vmm_v.at(j).second.hitsX.size() != 0)
        {
            vmm_hits_vec.clear();
            for (auto it = hits_vmm_v.at(j).second.hitsX.begin(); it != hits_vmm_v.at(j).second.hitsX.end(); ++it)
            {
                int strip = it->first * (1 - 8e-3) - 8.46 / 0.25;
                int pdo = it->second;
                if (pdo < 50)
                    continue;
                vmm_hits_vec.push_back(make_pair(strip, pdo));
            }
            if (vmm_hits_vec.size() == 0)
                continue;
            else
                num++;
        }
    }
    std::cout << "Num of good VMM hits: " << num << "\n";
    return 0;

    for (unsigned long j = 0; j < hits_vmm_v.size(); j++)
    {
        if (hits_vmm_v.at(j).second.sync && j >= 1485)
        {
            if (index == 0)
            {
                prevSyncBcid = hits_vmm_v.at(j).second.bcid;
                // std::cout << "VMM starrt t = " << hits_vmm_v.at(j).second.timeFull() << "\n";
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
                    if (!(diffDiff >= 2000 - 5 && diffDiff <= 2000 + 5) && !(diffDiff >= 4000 - 5 && diffDiff <= 4000 + 5))
                    {
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

        if (hits_vmm_v.at(j).second.hitsX.size() != 0)
        {
            vmm_hits_vec.clear();
            for (auto it = hits_vmm_v.at(j).second.hitsX.begin(); it != hits_vmm_v.at(j).second.hitsX.end(); ++it)
            {
                int strip = it->first * (1 - 8e-3) - 8.46 / 0.25;
                int pdo = it->second;
                if (pdo < 50)
                    continue;
                vmm_hits_vec.push_back(make_pair(strip, pdo));
            }
            if (vmm_hits_vec.size() == 0)
                continue;

            int diff_hit = 0;
            if (hits_vmm_v.at(j).second.bcid - prevSyncBcid >= 0)
            {
                diff_hit = hits_vmm_v.at(j).second.bcid - prevSyncBcid;
            }
            else
            {
                diff_hit = hits_vmm_v.at(j).second.bcid + 4096 - prevSyncBcid;
            }

            vmmHitTime = pulseTime + round(diff_hit * 25.0 / 1000.0);

            long long startT_apv = 0;
            long long startT_pulse_apv = 0;
            long long T_apv = 0;
            int nPeriodsAPV = 0;

            int prevSRS = -1;
            int prev_pulse_SRS = -1;
            long long pulser_T = 0;
            long long apv_hit_T = 0;
            bool hitMapped = false;

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
                        // std::cout << "APV starrt t = " << hits_apv_v.at(i).second.timeFull() << "\n";
                    }
                    else
                    {
                        nPeriodsAPV = round((T_apv - startT_pulse_apv) * 1.0 / 1e4);
                        prev_pulse_SRS = hits_apv_v.at(i).second.timeSrs;
                        pulser_T = nPeriodsAPV * 1e4;
                    }
                }

                if (prev_pulse_SRS != -1)
                {
                    apv_hit_T = pulser_T + apv_time_from_SRS(prev_pulse_SRS, hits_apv_v.at(i).second.timeSrs);
                    if (abs(apv_hit_T - vmmHitTime) > 1000 || j < 1485)
                        continue;

                    if (hits_apv_v.at(i).second.hitsPerLayer.at(0).size() != 0)
                    {
                        vector<pair<int, int>> apv_hits_vec;
                        vector<pair<int, int>> apv_hits_vec_l0;
                        vector<pair<int, int>> apv_hits_vec_l1;

                        for (auto &h : hits_apv_v.at(i).second.hitsPerLayer.at(0))
                        {
                            int strip = h.first;
                            int pdo = h.second;
                            if (strip > 118 && strip < 172)
                                apv_hits_vec_l0.push_back(make_pair(strip, pdo));
                        }
                        for (auto &h : hits_apv_v.at(i).second.hitsPerLayer.at(1))
                        {
                            int strip = h.first * (1 - 2.29e-3) - 2.412 / 0.25;
                            int pdo = h.second;
                            if (strip > 118 && strip < 172)
                                apv_hits_vec_l1.push_back(make_pair(strip, pdo));
                        }
                        for (auto &h : hits_apv_v.at(i).second.hitsPerLayer.at(2))
                        {
                            int strip = h.first * (1 - 8e-3) - 8.46 / 0.25;
                            int pdo = h.second;
                            if (strip > 118 && strip < 172)
                                apv_hits_vec.push_back(make_pair(strip, pdo));
                        }

                        if (apv_hits_vec.size() == 0)
                            continue;

                        hitMapped = false;

                        for (int k = 0; k < apv_hits_vec.size(); k++)
                        {
                            for (int l = 0; l < vmm_hits_vec.size(); l++)
                            {
                                if (hitMapped)
                                    continue;

                                if (vmm_hits_vec.size() != 0 && abs(apv_hits_vec.at(k).first - vmm_hits_vec.at(l).first) < 5)
                                    hitMapped = true;
                            }
                        }

                        if (hitMapped)
                        {
                            numOfMapped++;
                            std::cout << "VMM event: " << j << "\t Total of mapped " << numOfMapped << "\n";
                            if (vmmPrevHitT_mapped == 0)
                            {
                                vmmPrevHitT_mapped = vmmHitTime;
                            }
                            else
                            {
                                if (vmmHitTime - vmmPrevHitT_mapped == 0)
                                    continue;

                                std::cout << "VMM delta T from prev mapped hit: " << vmmHitTime - vmmPrevHitT_mapped << " [us]\n";
                                vmmPrevHitT_mapped = vmmHitTime;
                            }

                            if (apvPrevHitT_mapped == 0)
                            {
                                apvPrevHitT_mapped = apv_hit_T;
                            }
                            else
                            {
                                std::cout << "APV delta T from prev mapped hit: " << apv_hit_T - apvPrevHitT_mapped << " [us]\n";
                                apvPrevHitT_mapped = apv_hit_T;
                            }
                            // std::cout << "!!! Mapped hits: " << numOfMapped << " !!! \n";
                            out_txt << "------- APV Double ReadOut Period " << nPeriodsAPV << " -------- T = " << T_apv - startT_pulse_apv << " [us] \n";

                            for (int k = 0; k < apv_hits_vec.size(); k++)
                            {
                                out_txt << k << "\t Strip: " << apv_hits_vec.at(k).first << "\n";
                                out_txt << "  \t PDO: " << apv_hits_vec.at(k).second << "\n";
                            }
                            out_txt << "------- APV Layer 0 -------- \n";
                            for (int k = 0; k < apv_hits_vec_l0.size(); k++)
                            {
                                out_txt << k << "\t Strip: " << apv_hits_vec_l0.at(k).first << "\n";
                                out_txt << "  \t PDO: " << apv_hits_vec_l0.at(k).second << "\n";
                            }
                            out_txt << "------- APV Layer 1 -------- \n";
                            for (int k = 0; k < apv_hits_vec_l1.size(); k++)
                            {
                                out_txt << k << "\t Strip: " << apv_hits_vec_l1.at(k).first << "\n";
                                out_txt << "  \t PDO: " << apv_hits_vec_l1.at(k).second << "\n";
                            }

                            out_txt << "------- VMM Period " << nPeriods / 200 << "  (" << nPeriods % 200 << ") -------- T = " << pulseTime + diff_hit * 25 / 1000 << "[us] (dT = " << apv_hit_T - vmmHitTime << " [us]) \n";

                            for (int l = 0; l < vmm_hits_vec.size(); l++)
                            {
                                out_txt << "\t Strip: " << vmm_hits_vec.at(l).first << "\n";
                                out_txt << "\t PDO: " << vmm_hits_vec.at(l).second << "\n";
                            }
                        }
                    }
                }
            }

            if (!hitMapped)
            {
                n_of_unmapped_becauese_of_T++;

                std::cout << "\t VMM UNMAPPED HIT delta T from prev mapped hit: " << vmmHitTime - vmmPrevHitT_mapped << " [us]\n";
                std::cout << "\t N of unmapped hits " << n_of_unmapped_becauese_of_T << "\n";
            }
        }
    }
}
