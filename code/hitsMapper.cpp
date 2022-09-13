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

int vmmRemoveFirstNPuslers(vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> &hits_vmm_v, int nRemoved = 637){
    int firstSyncBcid = 0;
    int index = -nRemoved; // 650; 665

    unsigned long j = 0;
    for (j = 0; j < hits_vmm_v.size() && index < 0; j++)
    {
        if (hits_vmm_v.at(j).second.sync)
            index++;
    }
    hits_vmm_v.erase(hits_vmm_v.begin(), hits_vmm_v.begin() + j-1); // remove first j elements: from 0 to j-1
    if (index == 0) // check if arways?
    {
        firstSyncBcid = hits_vmm_v.at(0).second.bcid;
        hits_vmm_v.erase(hits_vmm_v.begin(), hits_vmm_v.begin()); // remove first elements
    }
    else if(index < 0)
    {
        cout << "Removed all VMM events due to incorrect index value" << endl;
    }
    return firstSyncBcid;
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
    out_APV.open("../out/APV_hits_maped_0832_423.txt");

    ofstream out_VMM;
    out_VMM.open("../out/VMM_hits_0832_423_after.txt");

    ofstream out_VMM_hits;
    out_VMM_hits.open("../out/VMM_hits_UNmaped_0832_423.txt");

    int numOfMapped = 0;

    auto out = TFile::Open("../out/mapped_0832_423.root", "recreate");

    auto mappedEventNums = new TTree("mappedEvents", "");
    long long eventNumAPV, eventNumVMM;
    mappedEventNums->Branch("apv", &eventNumAPV);
    mappedEventNums->Branch("vmm", &eventNumVMM);

    auto stripsVMM = make_shared<TH1F>("stripsVMM", "stripsVMM", 360, 0, 360);
    auto mappedHitsPdo = make_shared<TH1F>("mappedHitsPdo", "mappedHitsPdo", 2000, 0, 2000);
    auto mappedHitsPdo_apv = make_shared<TH1F>("mappedHitsPdo_apv", "mappedHitsPdo_apv", 2000, 0, 2000);
    auto hitsPdo = make_shared<TH2F>("hitsPdo", "hitsPdo", 360, 0, 360, 2000, 0, 2000);
    auto hitsPdo_apv = make_shared<TH2F>("hitsPdo_apv", "hitsPdo_apv", 360, 0, 360, 2000, 0, 2000);
    auto hitsPdoAll = make_shared<TH1F>("hitsPdoAll", "hitsPdoAll", 2000, 0, 2000);
    auto hitsPdoAll_apv = make_shared<TH1F>("hitsPdoAll_apv", "hitsPdoAll_apv", 2000, 0, 2000);
    auto layer2_vs_layer0 = make_shared<TH2F>("layer2_vs_layer0", "layer2_vs_layer0; layer2 strip №; #Delta strips", 200, 0, 200, 100, 0, 100);
    auto tbh_vmm = make_shared<TH1F>("tbh_vmm", "#Delta t between events with hits VMM in us", 2000, 0, 100000);
    auto tbh_vmm_mapped = make_shared<TH1F>("tbh_vmm_mapped", "#Delta t between events with mapped hits VMM in us", 2000, 0, 100000);
    auto tbh_vmm_umapped = make_shared<TH1F>("tbh_vmm_umapped", "#Delta t between events with unmapped hits VMM in us", 2000, 0, 100000);
    auto tbh_apv_mapped = make_shared<TH1F>("tbh_apv_mapped", "#Delta t between events with mapped hits APV in us", 2000, 0, 100000);
    auto tbh_apv_umapped = make_shared<TH1F>("tbh_apv_umapped", "#Delta t between events with unmapped hits APV in us", 2000, 0, 100000);
    auto tb_vmm_apv = make_shared<TH1F>("tb_vmm_apv", "#Delta t between events with mapped hits VMM and APV in us", 2000, 0, 100000);
    auto tb_u_vmm_apv = make_shared<TH1F>("tb_u_vmm_apv", "#Delta t between events with unmapped hits VMM and APV in us", 2000, 0, 100000);
    auto tbh_apv = make_shared<TH1F>("tbh_apv", "#Delta t between events with hits APV in us", 2000, 0, 20000);
    int a = -100;
    int a_hat = -5;
    int b = -100;
    long long prevT = 0;
    for (unsigned long j = 1485; j < hits_vmm_v.size(); j++)
    {
        if (hits_vmm_v.at(j).second.sync)
        {
            if (a < 0)
            {
                if (a == -100)
                {
                    std::cout << "VMM first pule DAQ t = " << hits_vmm_v.at(j).second.timeFull() << "\n \t PDO = " << hits_vmm_v.at(j).second.pdo << "\n";
                    std::cout << "VMM first pule J = " << j << "\n";
                    prevT = hits_vmm_v.at(j).second.timeFull();
                    a++;
                    continue;
                }
                // else
                // {
                //     std::cout << "VMM pulses delta T = " << hits_vmm_v.at(j).second.timeFull() - prevT << "\n \t PDO = " << hits_vmm_v.at(j).second.pdo << "\n";
                //     a++;
                //     prevT = hits_vmm_v.at(j).second.timeFull();
                // }
            }
            // else if (a_hat < 0)
            // {
            //     std::cout << "VMM F*ck'n delta T = " << hits_vmm_v.at(j).second.timeFull() - prevT << "\n \t PDO = " << hits_vmm_v.at(j).second.pdo << "\n";
            //     prevT = hits_vmm_v.at(j).second.timeFull();
            // }

            a_hat++;
        }

        if (hits_vmm_v.at(j).second.hitsX.size() != 0)
        {
            bool flag = false;
            for (auto it = hits_vmm_v.at(j).second.hitsX.begin(); it != hits_vmm_v.at(j).second.hitsX.end(); ++it)
            {
                int strip = it->first * (1 - 8e-3) - 8.46 / 0.25;
                int pdo = it->second;

                if (pdo < 50)
                    continue;

                stripsVMM->Fill(strip);
                hitsPdo->Fill(strip, pdo);
                hitsPdoAll->Fill(pdo);
                flag = true;
            }
            if (flag)
            {
                tbh_vmm->Fill(hits_vmm_v.at(j).second.timeFull() - prevT);
                prevT = hits_vmm_v.at(j).second.timeFull();
            }
        }
    }
    prevT = 0;

    for (unsigned long i = 0; i < hits_apv_v.size(); i++)
    {
        if (hits_apv_v.at(i).second.sync && b < 0)
        {
            if (b < 0)
            {
                if (b == -100)
                {
                    std::cout << "APV first pule DAQ t = " << hits_apv_v.at(i).second.timeFull() << "\n";
                    prevT = hits_apv_v.at(i).second.timeFull();
                    b++;
                    continue;
                }
                else
                {
                    std::cout << "APV pulses delta T = " << (hits_apv_v.at(i).second.timeFull() - prevT) / 1e3 << "\n";
                    b++;
                    prevT = hits_apv_v.at(i).second.timeFull();
                }
            }
        }
        if (hits_apv_v.at(i).second.hitsPerLayer.at(2).size() != 0)
        {
            bool flag = false;
            for (auto &h : hits_apv_v.at(i).second.hitsPerLayer.at(2))
            {
                int strip = h.first * (1 - 8e-3) - 8.46 / 0.25;
                int pdo = h.second;

                if (strip < 118 || strip > 172)
                    continue;

                hitsPdo_apv->Fill(strip, pdo);
                hitsPdoAll_apv->Fill(pdo);
                flag = true;
            }
            if (flag)
            {
                tbh_apv->Fill(hits_apv_v.at(i).second.timeFull() - prevT);
                prevT = hits_apv_v.at(i).second.timeFull();
            }
        }
    }
    int nHits = 0;
    for (unsigned long j = 0; j < hits_vmm_v.size(); j++)
    {

        if (hits_vmm_v.at(j).second.hitsX.size() != 0)
        {
            out_VMM << "------- VMM event " << j << "\n";

            for (auto it = hits_vmm_v.at(j).second.hitsX.begin(); it != hits_vmm_v.at(j).second.hitsX.end(); ++it)
            {
                int strip = it->first * (1 - 8e-3) - 8.46 / 0.25;
                int pdo = it->second;

                if (pdo < 50)
                    continue;

                nHits++;

                out_VMM << "\t Strip: " << strip << "\n";
                out_VMM << "\t PDO: " << pdo << "\n";
            }
        }
    }

    out_VMM << "\t Lost: " << nHits << "\n";
    out_VMM.close();

    // out->Write();
    // out->Close();

    // return 0;

    std::cout << "\t ---> TOTAL N of VMM events " << hits_vmm_v.size() << "\n";

    /* vmm: Remove first -index pulse syncs */
    int firstSyncBcid = vmmRemoveFirstNPuslers(hits_vmm_v);
    tuple<unsigned long, int, int, int, long long> beforeLastPulserParameters = {0, firstSyncBcid, 0, 0, 0};
    auto beforeLastPulserParametersCurrent = beforeLastPulserParameters;

    int mappedHitsVMM = 0;
    int UNmappedHitsVMM = 0;

    long long prevT_apv = 0;
    long long prevT_vmm = 0;

    int bad, prevSyncBcid, prevPrevSyncBcid, nPeriods;

    long long pulseTime;
    bool hitMapped;

    vector<pair<int, int>> vmm_hits_vec;
    vector<pair<int, int>> apv_hits_vec;
    vector<pair<int, int>> apv_hits_vec_l0;
    vector<pair<int, int>> apv_hits_vec_l1;

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
            // std::cout << "Period " << nPeriodsAPV << "--- is sync! N = " << hits_apv_v.at(i).second.hitsPerLayer.at(0).size() << "\n";
        }
        else
        {
            if (prev_pulse_SRS != -1)
            {
                // std::cout << "Not sync! N = " << hits_apv_v.at(i).second.hitsPerLayer.at(0).size() << "\n";
            }
        }

        if (prev_pulse_SRS != -1)
        {
            // std::cout << "\t Total:" << prev_pulse_SRS << "\t" << hits_apv_v.at(i).second.timeSrs << "\n";

            apv_hit_T = pulser_T + apv_time_from_SRS(prev_pulse_SRS, hits_apv_v.at(i).second.timeSrs);
            if (hits_apv_v.at(i).second.hitsPerLayer.at(0).size() != 0)
            {
                apv_hits_vec.clear();
                apv_hits_vec_l0.clear();
                apv_hits_vec_l1.clear();
                // std::cout << "---> with hits:" << prev_pulse_SRS << "\t" << hits_apv_v.at(i).second.timeSrs << "\n";
                // out_APV << "------- APV Period " << nPeriodsAPV << " -------- \n";
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
                

                vmm_hits_vec.clear();

                bad = 0;
                prevSyncBcid = get<1>(beforeLastPulserParameters);
                prevPrevSyncBcid = get<2>(beforeLastPulserParameters);
                nPeriods = get<3>(beforeLastPulserParameters);

                pulseTime = get<4>(beforeLastPulserParameters);
                hitMapped = false;
                beforeLastPulserParametersCurrent = beforeLastPulserParameters;

                for (unsigned long j = get<0>(beforeLastPulserParameters); j < hits_vmm_v.size(); j++)
                {
                    if (hits_vmm_v.at(j).second.sync)
                    {
                        beforeLastPulserParametersCurrent = {j, prevSyncBcid, prevPrevSyncBcid, nPeriods, pulseTime};
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

                    if (hits_vmm_v.at(j).second.hitsX.size() == 0 || nPeriods / 200 < nPeriodsAPV - 1)
                        continue;
                    else if (nPeriods / 200 > nPeriodsAPV + 1)
                        break;

                    int diff_hit = 0;
                    if (hits_vmm_v.at(j).second.bcid - prevSyncBcid >= 0)
                    {
                        diff_hit = hits_vmm_v.at(j).second.bcid - prevSyncBcid;
                    }
                    else
                    {
                        diff_hit = hits_vmm_v.at(j).second.bcid + 4096 - prevSyncBcid;
                    }
                    int dt_apv_vmm = T_apv - startT_pulse_apv - pulseTime - round(diff_hit * 25.0 / 1000.0);;
                    // int dt_apv_vmm = T_apv - startT_pulse_apv - pulseTime;

                    // std::cout << nPeriods / 200 << " \t " << hits_vmm_v.at(j).second.hitsX.size() << "\n";
                    // out_APV << "------- VMM Period " << nPeriods / 200 << "  (" << nPeriods % 200 << ") -------- dT = " << dt_apv_vmm << "\n";
                    if (abs(dt_apv_vmm) > 1000)
                        continue;

                    for (auto it = hits_vmm_v.at(j).second.hitsX.begin(); it != hits_vmm_v.at(j).second.hitsX.end(); ++it)
                    {
                        int strip = it->first * (1 - 8e-3) - 8.46 / 0.25;
                        int pdo = it->second;
                        if (pdo < 50)
                            continue;
                        vmm_hits_vec.push_back(make_pair(strip, pdo));
                        // out_APV << "\t Strip: " << it->first << "\n";
                        // out_APV << "\t PDO: " << it->second << "\n";
                    }
                    if (vmm_hits_vec.size() == 0)
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
                        // std::cout << "!!! Mapped hits: " << numOfMapped << " !!! \n";
                        if (T_apv - startT_pulse_apv - prevT_apv < 10e3)
                        {
                            out_APV << "------- APV Double ReadOut Period shiiiish " << nPeriodsAPV << " -------- T = " << T_apv - startT_pulse_apv << " [us] \n";
                        }
                        else
                        {
                            out_APV << "------- APV Double ReadOut Period " << nPeriodsAPV << " -------- T = " << T_apv - startT_pulse_apv << " [us] \n";
                        }
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
                            mappedHitsVMM++;
                        }
                        tbh_apv_mapped->Fill(T_apv - startT_pulse_apv - prevT_apv);
                        tbh_vmm_mapped->Fill(pulseTime + diff_hit * 25 / 1000 - prevT_vmm);
                        tb_vmm_apv->Fill(abs(dt_apv_vmm));


                        prevT_apv = T_apv - startT_pulse_apv;
                        prevT_vmm = pulseTime + diff_hit * 25 / 1000;

                    }
                    else
                    {
                        if (vmm_hits_vec.size() != 0)
                        {
                            out_VMM_hits << "------- VMM Period " << nPeriods / 200 << "  (" << nPeriods % 200 << ") -------- T = " << pulseTime + diff_hit * 25 / 1000 << "[us] (dT = " << dt_apv_vmm << " [us]) \n";

                            for (int l = 0; l < vmm_hits_vec.size(); l++)
                            {
                                out_VMM_hits << "\t Strip: " << vmm_hits_vec.at(l).first << "\n";
                                out_VMM_hits << "\t PDO: " << vmm_hits_vec.at(l).second << "\n";
                                UNmappedHitsVMM++;
                            }
                        }
                    }

                    if (hitMapped){
                        beforeLastPulserParameters = beforeLastPulserParametersCurrent;
                        eventNumAPV = hits_apv_v.at(i).first;
                        eventNumVMM = hits_vmm_v.at(j).first;
                        mappedEventNums->Fill();
                        break;
                    }
                }
            }
        }
    }
    std::cout << "N of MAPPED hits in VMM " << mappedHitsVMM << "\n";
    out_VMM_hits << "N of unMAPPED hits in VMM " << UNmappedHitsVMM << "\n";

    out_APV.close();
    out_VMM_hits.close();
    out->Write();
    out->Close();
}
