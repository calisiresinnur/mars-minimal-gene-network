"""
Ekstremofil (Salinibacter ruber, iMB631) için tekli gen silme analizi --
mars_gen_silme.py'nin (B. subtilis) analoğu.

DERS (mars_gen_silme.py'deki büyük düzeltmeden): solver tolerance baştan
1e-9'a sabitlendi -- düşük büyüme oranlarında sahte-feasible sonuç riskini
önlemek için.

Senaryolar: referans (kendi rahat koşulu) + iki Mars şiddet noktası (t=0.10
-> ~%10 büyüme, t=0.02 -> ~%2 büyüme -- B. subtilis'teki "hafif marj" ve
"sıkı marj"ın kabaca nispi analogları), 3 bakım çarpanında (×1.5/×2/×3).
"""

import os

import pandas as pd
from cobra.flux_analysis import single_gene_deletion

from extremofil_fba import bakim_reaksiyonunu_bul, mars_kisitlarini_uygula, modeli_yukle, referans_ortami_uygula

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SONUC_KLASORU = os.path.join(PROJE_KOKU, "results")
SOLVER_TOLERANCE = 1e-9
ESANSIYELLIK_ESIGI = 0.01
ISLEMCI_SAYISI = 4

SERT = dict(o2=0.02, glc=0.005, org=0.0005, h2o=0.05)
ILIMLI = dict(o2=3.15, glc=1.0, org=0.1, h2o=1000.0)


def kisit_degerleri(t):
    o2 = SERT["o2"] + t * (ILIMLI["o2"] - SERT["o2"])
    glc = SERT["glc"] + t * (ILIMLI["glc"] - SERT["glc"])
    org = SERT["org"] + t * (ILIMLI["org"] - SERT["org"])
    h2o = SERT["h2o"] + t * (ILIMLI["h2o"] - SERT["h2o"])
    return o2, glc, org, h2o


MARS_SENARYOLARI = []
for bakim in [1.5, 2.0, 3.0]:
    for t, etiket in [(0.10, "orta"), (0.02, "sert")]:
        o2, glc, org, h2o = kisit_degerleri(t)
        MARS_SENARYOLARI.append(dict(
            etiket=f"Mars_bakim_x{bakim}_{etiket}", o2=o2, glc=glc, org=org, h2o=h2o, bakim_carpani=bakim,
        ))


def senaryo_calistir(etiket, kisit_uygula):
    model = modeli_yukle()
    model.solver.configuration.tolerances.feasibility = SOLVER_TOLERANCE
    if kisit_uygula is not None:
        kisit_uygula(model)

    wt = model.optimize()
    print(f"{etiket}: WT büyüme = {wt.objective_value:.8f} (durum: {model.solver.status})")

    sonuc = single_gene_deletion(model, processes=ISLEMCI_SAYISI)
    sonuc = sonuc.reset_index(drop=True)
    sonuc["gen_id"] = sonuc["ids"].apply(lambda s: next(iter(s)) if s else None)
    sonuc["senaryo"] = etiket
    sonuc["wt_buyume"] = wt.objective_value
    sonuc["oran"] = sonuc["growth"] / wt.objective_value
    sonuc["esansiyel"] = sonuc["oran"] < ESANSIYELLIK_ESIGI
    return sonuc[["gen_id", "senaryo", "growth", "wt_buyume", "oran", "status", "esansiyel"]]


def mars_kisiti(senaryo):
    def uygula(model):
        atpm = bakim_reaksiyonunu_bul(model)
        mars_kisitlarini_uygula(
            model, atpm, o2_ub=senaryo["o2"], glc_ub=senaryo["glc"], organik_ub=senaryo["org"],
            h2o_ub=senaryo["h2o"], bakim_carpani=senaryo["bakim_carpani"], bakim_taban=3.15, sessiz=True,
        )
    return uygula


def main():
    os.makedirs(SONUC_KLASORU, exist_ok=True)

    # DİKKAT: None geçmek YETMEZ -- iMB631'in ham SBML varsayılan ortamı
    # tüm amino asit/vitaminleri sınırsız (1000) veriyor ve anlamsız bir
    # büyüme oranına yol açıyor (13.88/saat). referans_ortami_uygula()
    # AÇIKÇA çağrılmalı (bkz. extremofil_fba.py'deki kalibrasyon notu).
    tum_sonuclar = [senaryo_calistir("Referans", referans_ortami_uygula)]
    for s in MARS_SENARYOLARI:
        tum_sonuclar.append(senaryo_calistir(s["etiket"], mars_kisiti(s)))

    df = pd.concat(tum_sonuclar, ignore_index=True)
    csv_yolu = os.path.join(SONUC_KLASORU, "extremofil_gen_silme_sonuclari.csv")
    df.to_csv(csv_yolu, index=False)
    print(f"\nHam sonuçlar kaydedildi: {csv_yolu} ({len(df)} satır)")

    pivot = df.pivot(index="gen_id", columns="senaryo", values="esansiyel")
    ref_esansiyel = pivot["Referans"]
    mars_kolonlari = [c for c in pivot.columns if c != "Referans"]

    yeni_esansiyel = pivot[(~ref_esansiyel) & pivot[mars_kolonlari].any(axis=1)]
    dispanse = pivot[ref_esansiyel & (~pivot[mars_kolonlari]).any(axis=1)]

    print("\n--- Özet ---")
    for kolon in pivot.columns:
        print(f"{kolon:26s}: {int(pivot[kolon].sum()):4d} esansiyel gen / {len(pivot)}")
    print(f"\nMars'a özgü YENİ esansiyel gen adayı: {len(yeni_esansiyel)}")
    print(f"Mars'ta esansiyellikten çıkan gen: {len(dispanse)}")


if __name__ == "__main__":
    main()
