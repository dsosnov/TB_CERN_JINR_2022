#include "apv.C"
#include "evBuilder.C"
#include <iostream>
#include <fstream>


unsigned long apv_sync_index_finder(int whichSync, vector<pair<unsigned long, apv::doubleReadoutHits>> hits_apv_v)
{
    int index = 0;
    for (unsigned long i = 0; i < hits_apv_v.size(); i++)
    {
        if (hits_apv_v.at(i).second.sync)
        {
            index++;
            if (index == whichSync)
            {
                return i;
            }
        }
    }
    return 0;
} 


unsigned long vmm_sync_index_finder(int whichSync, vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> hits_vmm_v)
{
    int index = 0;
    for (unsigned long i = 0; i < hits_vmm_v.size(); i++)
    {
        if (hits_vmm_v.at(i).second.sync)
        {
            index++;
            if (index == whichSync)
            {
                return i;
            }
        }
    }
    return 0;
} 

int apv_window_finder(vector<pair<unsigned long, apv::doubleReadoutHits>> hits_apv_v)
{
    for (int j = 0; j < 1000; j++)
    {
        int start = apv_sync_index_finder(j, hits_apv_v);
        int stop = apv_sync_index_finder(j+1, hits_apv_v);
        if (hits_apv_v.at(stop).second.timeFull() - hits_apv_v.at(start).second.timeFull() < 1.1e4 && (hits_apv_v.at(start).second.hits.size() != 0 || hits_apv_v.at(start+1).second.hits.size() != 0))
        {
            return j;
        }
    }
    return -1;
    
}

long long apv_sync_time(int syncNum, vector<pair<unsigned long, apv::doubleReadoutHits>> hits_apv_v)
{
    int index = 0;
    long long prevSyncT = 0;
    long long tInMs = 0;

    for (unsigned long i = 0; i < hits_apv_v.size(); i++)
    {
        if (hits_apv_v.at(i).second.sync)
        {
            if (index == 0)
            {
                prevSyncT = hits_apv_v.at(i).second.timeFull();
                index++;
            }
            else 
            {
                index++;
                std::cout << round ((hits_apv_v.at(i).second.timeFull() - prevSyncT) / 1e4) << "\n";
                tInMs += 10 * round ((hits_apv_v.at(i).second.timeFull() - prevSyncT) / 1e4);
                prevSyncT = hits_apv_v.at(i).second.timeFull();
                if (index == syncNum)
                {
                    return tInMs;
                }
                    
            }
        }
    }
    return 0;
}


long long vmm_sync_time(int syncNum, vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> hits_vmm_v)
{
    int index = 0;
    int prevSyncBcid = 0;
    long long tInMus = 0;

    for (unsigned long i = 0; i < hits_vmm_v.size(); i++)
    {
        if (hits_vmm_v.at(i).second.sync)
        {
            if (index == 0)
            {
                prevSyncBcid = hits_vmm_v.at(i).second.bcid;
                index++;
            }
            else 
            {
                index++;

                int diff = 0;
                if (hits_vmm_v.at(i).second.bcid - prevSyncBcid > 0)
                {
                    diff = hits_vmm_v.at(i).second.bcid - prevSyncBcid;
                }
                else
                {
                    diff = hits_vmm_v.at(i).second.bcid + 4096 - prevSyncBcid;
                }
                std::cout << prevSyncBcid << "\t" << hits_vmm_v.at(i).second.bcid << "\t" << diff << "\n";
                tInMus += 50 * round(diff / 2000);
                prevSyncBcid = hits_vmm_v.at(i).second.bcid;
                if (index == syncNum)
                {
                    return tInMus;
                }
                    
            }
        }
    }
    return 0;
}

int vmm_window_finder(long long apvT, vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> hits_vmm_v)
{
    std::cout << vmm_sync_time(29200, hits_vmm_v) << std::endl;
    for (int i = 29200; i < 1e6; i++)
    {
        if (abs(apvT * 1e3 - vmm_sync_time(i, hits_vmm_v)) <= 10)
            return i;

    
    }
    return -1;
}

void apv_vmm_ana()
{
    pair<string, string> run_pair = {"run_0832_cut", "run423_cut"};

    auto out = TFile::Open("../out/apv_vmm_ana.root", "recreate");

    auto hSigDtAPV = make_shared<TH1F>("hSigDtAPV", "hSigDtAPV", 1e5, 0, 10 * 1e4);
    auto hSigDtVMM = make_shared<TH1F>("hSigDtVMM", "hSigDtVMM", 100, 0, 10 * 50);
    auto hSigDbcidVMM = make_shared<TH1F>("hSigDbcidVMM", "hSigDbcidVMM", 2050, 0, 4100);
    auto hBcidVMM = make_shared<TH1F>("hBcidVMM", "hBcidVMM", 2050, 0, 4100);
    auto hPeriodsAPV = make_shared<TH1F>("hPeriodsAPV", "hPeriodsAPV", 500, 0, 500);
    auto spillsAPV = make_shared<TH1F>("spillsAPV", "spillsAPV", 1e3, 0, 1e3);
    auto spillsVMM = make_shared<TH1F>("spillsVMM", "spillsVMM", 1e3, 0, 1e3);
    auto timeCheck = make_shared<TH2F>("timeCheck", "timeCheck; T_{DAQ}, us; T_{DAQ} - T_{PULSE}, us", 5e3, 0, 5e6, 1300, -200, 1100);
    auto timeCheckAPV = make_shared<TH2F>("timeCheckAPV", "timeCheckAPV; T_{DAQ}, us; T_{DAQ} - T_{PULSE}, us", 5e3, 0, 5e6, 1300, -200, 1100);
    auto VMM_time_check = make_shared<TH1F>("VMM_time_check", "VMM_time_check", 120e3, 0, 6e6);

    auto timeCheck_res = make_shared<TH2F>("timeCheck_res", "timeCheck_res; T_{DAQ}, us; T_{DAQ}^{hit_{i}} - T_{DAQ}^{hit_{i-1}}, us", 5e3, 0, 5e6, 600, -100, 500);
    auto timeCheckAPV_res = make_shared<TH2F>("timeCheckAPV_res", "timeCheckAPV_res; T_{DAQ}, us; T_{DAQ}^{hit_{i}} - T_{DAQ}^{hit_{i-1}}, us", 5e3, 0, 5e6, 400, 0.98e4, 1.02e4);

    auto timeCheck_res_ns = make_shared<TH2F>("timeCheck_res_ns", "timeCheck_res_ns; T_{BCID}, us; T_{BCID}^{hit_{i}} - T_{BCID}^{hit_{i-1}}, us", 5e3, 0, 5e6, 600, -100, 500);
    auto timeCheckAPV_res_ns = make_shared<TH2F>("timeCheckAPV_res_ns", "timeCheckAPV_res_ns; T_{SRS}, us; T_{SRS}^{hit_{i}} - T_{SRS}^{hit_{i-1}}, us", 5e3, 0, 5e6, 400, 0.98e4, 1.02e4);

    unsigned long long from = 0, to = 0;

    auto apvan = new apv(run_pair.second);
    apvan->useSyncSignal();
    // auto hits_apv = apvan->GetCentralHits(from, to);
    auto hits_apv = apvan->GetCentralHits2ROnly(from, to);
    auto vmman = new evBuilder(run_pair.first, "g1_p25_s100-0&60", "map-20220605");
    vmman->useSyncSignal();
    auto hits_vmm = vmman->GetCentralHits(from, to);

    for (auto &h : hits_vmm)
        h.second.time += 325;

    vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> hits_vmm_v;
    hits_vmm_v.assign(hits_vmm.begin(), hits_vmm.end());
    // vector<pair<unsigned long, analysisGeneral::mm2CenterHitParameters>> hits_apv_v;
    vector<pair<unsigned long, apv::doubleReadoutHits>> hits_apv_v;
    hits_apv_v.assign(hits_apv.begin(), hits_apv.end());

    for (unsigned long i = 0; i < hits_apv_v.size(); i++)
    {
        if (!hits_apv_v.at(i).second.sync)
        {
            spillsAPV->Fill(hits_apv_v.at(i).second.timeSec - 1.654545e9);

        }
        else
            continue;
    }

    for (unsigned long i = 0; i < hits_vmm_v.size(); i++)
    {
        if (hits_vmm_v.at(i).second.trigger)
        {
            spillsVMM->Fill(hits_vmm_v.at(i).second.timeSec - 1.654545e9);

        }
        else
            continue;
    }

    long long prevSync = -1;

    for (int i = 0; i < hits_apv_v.size(); i++)
    {
        if (hits_apv_v.at(i).second.sync)
        {
            if (prevSync == -1)
            {
                prevSync = hits_apv_v.at(i).second.timeFull(); 
            }
            else
            {
                hSigDtAPV->Fill(hits_apv_v.at(i).second.timeFull() - prevSync);
                prevSync = hits_apv_v.at(i).second.timeFull();
            }
        }
        else
            continue;
    }
    int index = -5;
    int bad = 0;
    int prevSyncBcid = 0;
    int prevPrevSyncBcid = 0;
    int nPeriods = 0;
    long long startT_vmm = 0;
    long long stopT_vmm = 0;

    double pulseTime = 0;
    double daqTime = 0;
    double prevDaqTime = 0;
    double prevDaqTime_res = 0;

    long long bcidT = 0;

    std::cout << "start ------------- \n";

    ofstream out_VMM;
    out_VMM.open ("../out/VMM_hits.txt");
    for (unsigned long i = 0; i < hits_vmm_v.size(); i++)
    {
        // if (hits_vmm_v.at(i).second.timeSec - 1.654545e9 < 164 || hits_vmm_v.at(i).second.timeSec - 1.654545e9 > 169)
        // if (hits_vmm_v.at(i).second.timeSec - 1.654545e9 < 160 || hits_vmm_v.at(i).second.timeSec - 1.654545e9 > 163)
            // continue;
        
        if (hits_vmm_v.at(i).second.hitsX.size() != 0 && daqTime / 1e3 > 33)
        {   
            out_VMM << "------- VMM Period " << nPeriods / 200 - 3 << "  (" << nPeriods % 200 << ") -------- \n";
            for(std::map<unsigned int,unsigned int>::iterator it = hits_vmm_v.at(i).second.hitsX.begin(); it != hits_vmm_v.at(i).second.hitsX.end(); ++it){
                out_VMM << "\t Strip: " << it->first << "\n";
                out_VMM << "\t PDO: " << it->second << "\n";
            }
        }

        if (hits_vmm_v.at(i).second.sync)
        {
            if (index < 0)
            {
                index++;
                continue;
            }
            if (index == 0)
            {
                prevSyncBcid = hits_vmm_v.at(i).second.bcid;
                startT_vmm = hits_vmm_v.at(i).second.timeFull();
                index++;
            }
            else 
            {

                
                daqTime = hits_vmm_v.at(i).second.timeFull() - startT_vmm;
                VMM_time_check->Fill(daqTime - prevDaqTime);
                prevDaqTime = daqTime;
                index++;

                int diff = 0;
                if (hits_vmm_v.at(i).second.bcid - prevSyncBcid > 0)
                {
                    diff = hits_vmm_v.at(i).second.bcid - prevSyncBcid;
                    // std::cout << i << "\t" << diff << "\n";
                }
                else
                {
                    diff = hits_vmm_v.at(i).second.bcid + 4096 - prevSyncBcid;
                    // std::cout << i << "\t" << diff << "\n";
                }
                int diffDiff = 0; 
                
                if (hits_vmm_v.at(i).second.bcid - prevPrevSyncBcid > 0)
                {
                    diffDiff = hits_vmm_v.at(i).second.bcid - prevPrevSyncBcid;
                }
                else
                {
                    diffDiff = hits_vmm_v.at(i).second.bcid + 4096 - prevPrevSyncBcid;
                }

                if (daqTime < 200e3 && diff != 2000)
                { 
                    std::cout << daqTime/1e3 << "\t" << diff << "\n";
                }

                

                if(!(diff >= 2000 - 5 && diff <= 2000 + 5) && !(diff >= 4000 - 5 && diff <= 4000 + 5))
                {
                    if (diff == 1904)
                    {
                        // std::cout << i << " \t ----> " << diff << "\t DAQT " << daqTime << "\n";
                        diff = 6000;
                    }
                    else if (diff == 3904)
                    {
                        // std::cout << i << " \t ----> " << diff << "\t DAQT " << daqTime << "\n";
                        diff = 8000;
                    }
                    else if (!(diffDiff >= 2000 - 5 && diffDiff <= 2000 + 5) && !(diffDiff >= 4000 - 5 && diffDiff <= 4000 + 5))
                    {
                        // std::cout << i << "\t" << prevSyncBcid << "\t" << hits_vmm_v.at(i).second.bcid << "\t" << diff << "\n";
                        bad++;
                        prevPrevSyncBcid = prevSyncBcid;
                        prevSyncBcid = hits_vmm_v.at(i).second.bcid;
                        continue;
                    }
                    else
                    {
                        diff = diffDiff;
                    }
                }
                bcidT += diff * 25 / 1000;
                timeCheck_res_ns->Fill(bcidT, diff * 25 / 1000);
                
                // if (nPeriods == 200 * 10)
                //     break;

                prevPrevSyncBcid = prevSyncBcid;
                prevSyncBcid = hits_vmm_v.at(i).second.bcid;
                nPeriods += round (diff * 1.0 / 2000.0);
                pulseTime = nPeriods * 50;
                stopT_vmm = hits_vmm_v.at(i).second.timeFull();

                

                timeCheck->Fill(daqTime, daqTime - bcidT);
                // out_VMM << (daqTime - prevDaqTime_res) / 1e3 << "\t" << daqTime / 1e6 << "\n";
                timeCheck_res->Fill(daqTime, daqTime - prevDaqTime_res);
                prevDaqTime_res = daqTime;
         
            }
        }
        else 
        {
            // spillsVMM->Fill(hits_vmm_v.at(i).second.timeFull() / 1e6);
        }
    }
    out_VMM.close();
    std::cout << "BAD ------------- " << bad << "\n";
    std::cout << "VMM periods ------------- " << nPeriods / 200 << "\n";
    std::cout << "VMM strart " << startT_vmm << "\t stop " << stopT_vmm << "\t dif " << (stopT_vmm - startT_vmm) / 1e6 << "\t PulseTIME " << pulseTime / 1e6 << "\n";

    int indexAPV = 0;
    long long prevSyncT = 0;
    long long startT_apv = 0;
    long long startT_pulse_apv = 0;
    long long stopT_apv = 0;
    int nPeriodsAPV = 0;
    int srsT_period = 16777215;
    long long int srs_tempT = 0;

    long long startT_DAQ = 0;
    long long stopT_DAQ = 0;
    long long prevT_DAQ = 0;
    int prevSRS = 0;

    ofstream out_APV;
    out_APV.open ("../out/APV_hits.txt");

    for (unsigned long i = 0; i < hits_apv_v.size(); i++)
    {
        // if (hits_apv_v.at(i).second.timeSec - 1.654545e9 < 164 || hits_apv_v.at(i).second.timeSec - 1.654545e9 > 169)
        // if (hits_apv_v.at(i).second.timeSec - 1.654545e9 < 160 || hits_apv_v.at(i).second.timeSec - 1.654545e9 > 163)
            // continue;

        if (hits_apv_v.at(i).second.hits.size() != 0)    
        { 
            out_APV << "------- APV Period " << round ((srs_tempT - startT_pulse_apv) * 1.0 / 1e7) << " -------- \n";
            for (int j = 0; j < hits_apv_v.at(i).second.hits.size(); j++)
            {
                out_APV << j << "\t Strip: " << hits_apv_v.at(i).second.hits.at(j).strip << "\n";
                out_APV << "  \t PDO: " << hits_apv_v.at(i).second.hits.at(j).max_q << "\n";
            }
        }
        
        if (prevSyncT == 0)
        {
            prevSyncT = hits_apv_v.at(i).second.timeSrs;
            startT_apv = 25 * hits_apv_v.at(i).second.timeSrs;
        }
        else
        {
            if (prevSyncT > hits_apv_v.at(i).second.timeSrs)
            {
                srs_tempT += 25 * (hits_apv_v.at(i).second.timeSrs + srsT_period - prevSyncT);
            }
            else
            {
                srs_tempT += 25 * (hits_apv_v.at(i).second.timeSrs - prevSyncT);
            }
            prevSyncT = hits_apv_v.at(i).second.timeSrs;
        }   

        if (hits_apv_v.at(i).second.sync)
        {
            if (indexAPV == 0)
            {
                prevSRS = hits_apv_v.at(i).second.timeSrs;
                indexAPV++;
                startT_pulse_apv = srs_tempT;
                startT_DAQ = hits_apv_v.at(i).second.timeFull();
                stopT_DAQ = hits_apv_v.at(i).second.timeFull();
            }
            else 
            {
                prevT_DAQ = (stopT_DAQ - startT_DAQ);
                stopT_DAQ = hits_apv_v.at(i).second.timeFull();
                timeCheckAPV_res->Fill(stopT_DAQ - startT_DAQ, stopT_DAQ - startT_DAQ - prevT_DAQ);
                timeCheckAPV_res_ns->Fill(srs_tempT / 1e3, srs_tempT / 1e3 - stopT_apv / 1e3);
                prevSRS = hits_apv_v.at(i).second.timeSrs;
                
                indexAPV++;
                nPeriodsAPV = round ((srs_tempT - startT_pulse_apv) * 1.0 / 1e7);

                stopT_apv = srs_tempT - startT_pulse_apv;
                timeCheckAPV->Fill(stopT_DAQ - startT_DAQ, stopT_DAQ - startT_DAQ - srs_tempT / 1e3);
                // if (nPeriodsAPV == 10)
                //     break;
            }
        }
        else
        {
            // spillsAPV->Fill(hits_apv_v.at(i).second.timeFull() / 1e6);
        }
    }
    out_APV.close();
    std::cout << "APV periods ------------- " << nPeriodsAPV << "\n";
    std::cout << "APV strart " << startT_pulse_apv << "\t stop " << stopT_apv << "\t dif " << stopT_apv / 1e9 << "\t PulseTIME " << nPeriodsAPV * 0.01 << "\n";
    std::cout << "APV strart DAQ " << startT_DAQ << "\t stop DAQ " << stopT_DAQ << "\t dif " << (stopT_DAQ - startT_DAQ) / 1e6 << "\t PulseTIME " << nPeriodsAPV * 0.01 << "\n";
    out->Write();
    out->Close();

    auto canvas = new TCanvas("canvas", "canvas", 1500, 900);
    canvas->cd();
    gStyle->SetOptStat(0);
    timeCheck->SetLineColor(kBlack);
    timeCheck->Draw();
    timeCheckAPV->SetLineColor(kMagenta);
    timeCheckAPV->SetMarkerColor(kMagenta);
    timeCheckAPV->SetMarkerStyle(21);
    timeCheckAPV->SetMarkerSize(1.1);
    timeCheckAPV->Draw("SAME");
    auto legend = new TLegend(0.75, 0.3, 0.9, 0.15, "");
    legend->AddEntry("timeCheck", "VMM", "l");
    legend->AddEntry("timeCheckAPV", "APV", "l");
    legend->Draw();
    canvas->SaveAs("../out/time_check.pdf");

    return;

    // for (int i = 1; i < hits_vmm_v.size(); i++)
    // {
    //     // if (hits_vmm_v.at(i).second.sync && hits_vmm_v.at(i-1).second.sync)
    //     if (true)
    //     {
    //         // hBcidVMM->Fill(hits_vmm_v.at(i).second.bcid);
    //         hSigDtVMM->Fill(hits_vmm_v.at(i).second.timeFull() - hits_vmm_v.at(i - 1).second.timeFull());
    //         // hSigDbcidVMM->Fill(hits_vmm_v.at(i).second.bcid - hits_vmm_v.at(i - 1).second.bcid);
    //     }
    // }

    int apv_window_start = apv_window_finder(hits_apv_v);
    long long apv_t_start = apv_sync_time(apv_window_start, hits_apv_v);
    std::cout << "APV num of period: " << apv_window_start  << "\t  T_start[ms]: "<< apv_t_start << "\n";
    unsigned long apv_start_i = apv_sync_index_finder(apv_window_start, hits_apv_v);
    unsigned long apv_stop_i = apv_sync_index_finder(apv_window_start + 1, hits_apv_v);
    int vmm_window_start = vmm_window_finder(apv_t_start, hits_vmm_v);
    unsigned long vmm_start_i = vmm_sync_index_finder(vmm_window_start, hits_vmm_v);
    unsigned long vmm_stop_i = vmm_sync_index_finder(vmm_window_start + 200, hits_vmm_v);
    long long vmmSyncT = vmm_sync_time(vmm_window_start, hits_vmm_v);
    std::cout << "VMM num of period: " << vmm_window_start  << "\t  T_start[mus]: "<< vmmSyncT << "\n";
    
    // for(int i = 0; i < hits_apv_v.at(apv_start_i).second.hitsSync.size(); i++)
    // {
    //     std::cout << i << "\t sync " << hits_apv_v.at(apv_start_i).second.hitsSync.at(i).strip << "\t" << apv_t_start << "\n";
    // }
    
    std::cout << "-------------------------- \n APV: \n";

    for (unsigned long i = apv_start_i; i < apv_stop_i; i++)
    {
        if (hits_apv_v.at(i).second.hits.size() == 0 && hits_apv_v.at(i+1).second.hits.size() == 0)
            continue;

        for (int j = 0; j < hits_apv_v.at(i).second.hits.size(); j++)
        {
            std::cout << j << "\t Strip: " << hits_apv_v.at(i).second.hits.at(j).strip << std::endl;
            std::cout << "  \t PDO: " << hits_apv_v.at(i).second.hits.at(j).max_q << std::endl;
        }
    }
    std::cout << "-------------------------- \n VMM \n";
    int numOfSyncs = 0;
    long long timeFromStartSync = 0;
    // for (unsigned long i = vmm_start_i; i < vmm_stop_i; i++)
    // {
    //     long long vmmT = vmm_sync_time(i, hits_vmm_v);
    //     timeFromStartSync += vmmT - vmmSyncT;
    //     // if (timeFromStartSync > 1e4)
    //     //     break;
    //     std::cout << i << "\t" << timeFromStartSync << "\t" << vmmT << "\n";
    //     if (hits_vmm_v.at(i).second.sync && hits_vmm_v.at(i).second.hitsX.size() == 0)
    //     {
    //         std::cout << i << "\t sync only \t" << timeFromStartSync << "\t" << vmmT << "\n";
    //         numOfSyncs++;
    //     }
            
    //     else if (hits_vmm_v.at(i).second.sync && hits_vmm_v.at(i).second.hitsX.size() != 0)
    //     {
    //         std::cout << i << "\t sync with hits \t" << timeFromStartSync << "\t" << vmmT << "\n";
    //         for (auto const& imap: hits_vmm_v.at(i).second.hitsX)
    //         {
    //             std::cout << "\t Strip: " << imap.first << std::endl;
    //             std::cout << "\t PDO: " << imap.second << std::endl;
    //         }
    //         numOfSyncs++;
    //     }

    //     else if (!hits_vmm_v.at(i).second.sync && hits_vmm_v.at(i).second.hitsX.size() != 0)
    //     {
    //         std::cout << i << "\t just hits \t" << timeFromStartSync << "\t" << vmmT << "\n";
    //         for (auto const& imap: hits_vmm_v.at(i).second.hitsX)
    //         {
    //             std::cout << "\t Strip: " << imap.first << std::endl;
    //             std::cout << "\t PDO: " << imap.second << std::endl;
    //         }
    //     }
    //     else
    //     {
    //         continue;
    //     }
    // }

    for (unsigned long i = vmm_start_i; i < vmm_start_i + 1000; i++)
    {
        long long vmmT = vmm_sync_time(vmm_window_start + numOfSyncs, hits_vmm_v);
        if (hits_vmm_v.at(i).second.sync)
        {
            timeFromStartSync = vmmT - vmmSyncT;
            if (timeFromStartSync > 1e4)
                break;
            std::cout << i << "\t sync only \t" << timeFromStartSync << "\t" << vmmT << "\n";
            numOfSyncs++;
            if (hits_vmm_v.at(i).second.hitsX.size() != 0)
            {
                std::cout << i << "\t sync with hits \t" << timeFromStartSync << "\t" << vmmT << "\n";
                for (auto const& imap: hits_vmm_v.at(i).second.hitsX)
                {
                    std::cout << "\t Strip: " << imap.first << std::endl;
                    std::cout << "\t PDO: " << imap.second << std::endl;
                }
            }

        } 
        else if (hits_vmm_v.at(i).second.hitsX.size() != 0)
        {
            std::cout << i << "\t just hits \t" << timeFromStartSync << "\t" << vmmT << "\n";
            for (auto const &imap : hits_vmm_v.at(i).second.hitsX)
            {
                std::cout << "\t Strip: " << imap.first << std::endl;
                std::cout << "\t PDO: " << imap.second << std::endl;
            }
        }
        else 
            continue;
    }
    


    std::cout << "Num sync signals in VMM \t" << numOfSyncs << "\n";

    out->Write();
    out->Close();
}
