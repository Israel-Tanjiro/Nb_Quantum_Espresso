// ================================================================
//  Dos_Paper.C — Phonon Branch-Resolved DOS for Nb (BCC)
//  Reads per-branch files from branch_dos.py output
//  Usage: root -l 'Dos_Paper.C'
// ================================================================

void Dos_Paper() {

    // ============================================================
    //  USER SETTINGS
    // ============================================================
    const Int_t    N_POINTS   = 2000;       // must match N_W in branch_dos.py
    const Double_t X_MIN      = 0.0;        // THz
    const Double_t X_MAX      = 7.6;        // THz
    const Double_t CONV_MEV   = 4.14;       // 1 THz = 4.14 meV
    const Double_t FIT_MIN    = 0.8;        // THz — cubic fit start
    const Double_t FIT_MAX    = 2.0;        // THz — cubic fit end
    const Double_t EVAL_AT    = 1.0;        // THz — evaluate fit here

    TString FILE_TA1   = "DOS_Full_q_real_TA1.txt";
    TString FILE_TA2   = "DOS_Full_q_real_TA2.txt";
    TString FILE_LA    = "DOS_Full_q_real_LA.txt";
    TString FILE_TOTAL = "DOS_Full_q_real_Total.txt";

    // ============================================================
    //  READ DATA — one file per branch
    // ============================================================
    auto ReadBranch = [&](TString filename, Int_t n) {
        Float_t *freq = new Float_t[n];
        Float_t *dos  = new Float_t[n];
        ifstream f(filename);
        for (Int_t i = 0; i < n; i++) {
            f >> freq[i] >> dos[i];
            if (!f.good()) { cout << "ERROR reading " << filename << endl; break; }
        }
        f.close();
        return make_pair(freq, dos);
    };

    auto [freq,    dos_total] = ReadBranch(FILE_TOTAL, N_POINTS);
    auto [freq_1,  dos_1]     = ReadBranch(FILE_TA1,   N_POINTS);
    auto [freq_2,  dos_2]     = ReadBranch(FILE_TA2,   N_POINTS);
    auto [freq_3,  dos_3]     = ReadBranch(FILE_LA,    N_POINTS);

    // ============================================================
    //  CANVAS 1 — Full DOS
    // ============================================================
    TCanvas *c1 = new TCanvas("c1", "Phonon DOS — Nb BCC", 900, 700);
    c1->SetLeftMargin(0.12);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.12);

    TMultiGraph *mg = new TMultiGraph();

    TGraph *gTotal = new TGraph(N_POINTS, freq,   dos_total);
    TGraph *gTA1   = new TGraph(N_POINTS, freq_1, dos_1);
    TGraph *gTA2   = new TGraph(N_POINTS, freq_2, dos_2);
    TGraph *gLA    = new TGraph(N_POINTS, freq_3, dos_3);

    // Style — Total
    gTotal->SetLineColor(kBlack);
    gTotal->SetLineStyle(7);
    gTotal->SetLineWidth(3);
    gTotal->SetTitle("Total");

    // Style — TA1
    gTA1->SetLineColor(kGreen+1);
    gTA1->SetLineStyle(1);
    gTA1->SetLineWidth(2);
    gTA1->SetTitle("TA_{1}");

    // Style — TA2
    gTA2->SetLineColor(kRed);
    gTA2->SetLineStyle(1);
    gTA2->SetLineWidth(2);
    gTA2->SetTitle("TA_{2}");

    // Style — LA
    gLA->SetLineColor(kBlue);
    gLA->SetLineStyle(1);
    gLA->SetLineWidth(2);
    gLA->SetTitle("LA");

    mg->Add(gTA1);
    mg->Add(gTA2);
    mg->Add(gLA);
    mg->Add(gTotal);
    mg->Draw("AC");

    mg->GetXaxis()->SetTitle("Energy (THz)");
    mg->GetYaxis()->SetTitle("Density of States [a.u]");
    mg->GetXaxis()->SetTitleSize(0.05);
    mg->GetYaxis()->SetTitleSize(0.05);
    mg->GetXaxis()->SetLabelSize(0.05);
    mg->GetYaxis()->SetLabelSize(0.05);
    mg->GetXaxis()->CenterTitle(true);
    mg->GetYaxis()->CenterTitle(true);
    mg->GetXaxis()->SetLimits(X_MIN, X_MAX);

    // Legend
    TLegend *leg = new TLegend(0.65, 0.65, 0.92, 0.88);
    leg->AddEntry(gTotal, "Total DOS",  "l");
    leg->AddEntry(gTA1,   "TS", "l");
    leg->AddEntry(gTA2,   "TF", "l");
    leg->AddEntry(gLA,    "LA",     "l");
    leg->SetTextSize(0.04);
    leg->Draw();

    // meV axis on top
    Double_t ymax   = gPad->GetUymax();
    Double_t xmax_mev = X_MAX * CONV_MEV;
    TF1 *fconv = new TF1("fconv", Form("%f*x", CONV_MEV), 0, xmax_mev);
    TGaxis *axisMeV = new TGaxis(X_MIN, ymax, X_MAX, ymax,
                                  "fconv", 510, "-");
    axisMeV->SetTitle("Energy (meV)");
    axisMeV->CenterTitle(true);
    axisMeV->SetLabelSize(0.05);
    axisMeV->SetTitleSize(0.05);
    axisMeV->SetLabelFont(42);
    axisMeV->SetTitleFont(42);
    axisMeV->SetTitleOffset(1.2);
    axisMeV->Draw();
    leg->Draw();

    // Max frequency annotation
    Double_t max_freq = *max_element(freq, freq + N_POINTS);
    TPaveText *pt = new TPaveText(max_freq - 0.3, ymax * 0.1,
                                   max_freq + 0.2, ymax * 0.25, "NB");
    pt->AddText("#omega_{A}");
    pt->SetFillColor(0);
    pt->SetBorderSize(0);
    pt->Draw();

    TArrow *ar = new TArrow(max_freq, ymax*0.05,
                             max_freq, ymax*0.3, 0.02, "<|");
    ar->SetAngle(60);
    ar->SetLineWidth(2);
    ar->SetArrowSize(0.02);
    ar->Draw();

    c1->SaveAs("Nb_DOS_full.pdf");

    // ============================================================
    //  CANVAS 2 — Fractional DOS (0 to 2 THz) + cubic fit
    // ============================================================

    // Find index where freq < 2.0 THz
    Int_t j = 0;
    while (j < N_POINTS && freq[j] < 2.0) j++;

    Float_t *FreqF  = new Float_t[j];
    Float_t *Frac1  = new Float_t[j];
    Float_t *Frac2  = new Float_t[j];
    Float_t *Frac3  = new Float_t[j];

    for (Int_t i = 0; i < j; i++) {
        FreqF[i] = freq[i];
        Double_t tot = dos_1[i] + dos_2[i] + dos_3[i];
        if (tot > 1e-10) {
            Frac1[i] = dos_1[i] / tot;
            Frac2[i] = dos_2[i] / tot;
            Frac3[i] = dos_3[i] / tot;
        } else {
            Frac1[i] = Frac2[i] = Frac3[i] = 0;
        }
    }

    TCanvas *c2 = new TCanvas("c2", "Fractional DOS — Nb BCC", 900, 700);
    c2->SetLeftMargin(0.12);
    c2->SetRightMargin(0.05);

    TMultiGraph *mgF = new TMultiGraph();

    TGraph *gF1 = new TGraph(j, FreqF, Frac1);
    TGraph *gF2 = new TGraph(j, FreqF, Frac2);
    TGraph *gF3 = new TGraph(j, FreqF, Frac3);

    gF1->SetLineColor(kGreen+1); gF1->SetLineWidth(3);
    gF1->SetLineStyle(7);        gF1->SetTitle("TA_{1}");

    gF2->SetLineColor(kCyan+1);  gF2->SetLineWidth(3);
    gF2->SetLineStyle(7);        gF2->SetTitle("TA_{2}");

    gF3->SetLineColor(kRed);     gF3->SetLineWidth(3);
    gF3->SetLineStyle(7);        gF3->SetTitle("LA");

    // Cubic fits (pol3) — matches ROOT convention
    TF1 *fTA1 = new TF1("fTA1", "pol3(0)", FIT_MIN, FIT_MAX);
    TF1 *fTA2 = new TF1("fTA2", "pol3(0)", FIT_MIN, FIT_MAX);
    TF1 *fLA  = new TF1("fLA",  "pol3(0)", FIT_MIN, FIT_MAX);

    fTA1->SetLineColor(kGreen+1); fTA1->SetLineWidth(2);
    fTA2->SetLineColor(kCyan+1);  fTA2->SetLineWidth(2);
    fLA->SetLineColor(kRed);      fLA->SetLineWidth(2);

    gF1->Fit("fTA1", "Q", "", FIT_MIN, FIT_MAX);
    gF2->Fit("fTA2", "Q", "", FIT_MIN, FIT_MAX);
    gF3->Fit("fLA",  "Q", "", FIT_MIN, FIT_MAX);

    cout << "\n  ================================================" << endl;
    cout << "  Fractional DOS at " << EVAL_AT << " THz (cubic fit):" << endl;
    cout << "  TA1 = " << fTA1->Eval(EVAL_AT) << endl;
    cout << "  TA2 = " << fTA2->Eval(EVAL_AT) << endl;
    cout << "  LA  = " << fLA->Eval(EVAL_AT)  << endl;
    cout << "  Sum = " << fTA1->Eval(EVAL_AT) +
                          fTA2->Eval(EVAL_AT) +
                          fLA->Eval(EVAL_AT)  << endl;
    cout << "  ================================================\n" << endl;

    mgF->Add(gF1);
    mgF->Add(gF2);
    mgF->Add(gF3);
    mgF->Draw("AC");

    mgF->GetXaxis()->SetTitle("Energy (THz)");
    mgF->GetYaxis()->SetTitle("Fractional DOS");
    mgF->GetXaxis()->SetTitleSize(0.05);
    mgF->GetYaxis()->SetTitleSize(0.05);
    mgF->GetXaxis()->SetLabelSize(0.05);
    mgF->GetYaxis()->SetLabelSize(0.05);
    mgF->GetXaxis()->CenterTitle(true);
    mgF->GetYaxis()->CenterTitle(true);
    mgF->GetXaxis()->SetLimits(0.0, 2.0);
    mgF->GetYaxis()->SetRangeUser(0.0, 1.05);

    // Reference line at 1/3
    TLine *lRef = new TLine(0.0, 1.0/3.0, 2.0, 1.0/3.0);
    lRef->SetLineStyle(2);
    lRef->SetLineColor(kGray+1);
    lRef->SetLineWidth(1);
    lRef->Draw();

    // Mark evaluation point
    TLine *lEval = new TLine(EVAL_AT, 0.0, EVAL_AT, 1.05);
    lEval->SetLineStyle(3);
    lEval->SetLineColor(kBlack);
    lEval->SetLineWidth(2);
    lEval->Draw();

    TLegend *legF = new TLegend(0.65, 0.65, 0.92, 0.88);
    legF->AddEntry(gF1,  "TS", "l");
    legF->AddEntry(gF2,  "TF", "l");
    legF->AddEntry(gF3,  "LA",     "l");
    legF->AddEntry(lRef, "1/3",    "l");
    legF->SetTextSize(0.04);
    legF->Draw();

    c2->SaveAs("Nb_DOS_fractional.pdf");

    // Clean up
    delete[] freq;    delete[] dos_total;
    delete[] freq_1;  delete[] dos_1;
    delete[] freq_2;  delete[] dos_2;
    delete[] freq_3;  delete[] dos_3;
    delete[] FreqF;
    delete[] Frac1;   delete[] Frac2;   delete[] Frac3;
}
