"""
Ekstremofil (Salinibacter ruber, iMB631) için Mars kısıt şiddetine duyarlılık
analizi -- mars_duyarlilik.py'nin (B. subtilis) doğrudan analoğu.

BAŞTAN SIKI SOLVER TOLERANCE (kritik ders, mars_gen_silme.py'deki büyük
düzeltmeden öğrenildi): Çok düşük büyüme oranlarında biyokütle denklemindeki
kofaktör akışları solver'ın varsayılan toleransına (1e-7) yakınsayabilir ve
sahte-feasible sonuçlar üretebilir. Bu script baştan itibaren 1e-9 kullanıyor.

Şiddet ekseni ankorları burada B. subtilis'ten TAMAMEN FARKLI mutlak sayılar
-- iki modelin ölçeği/birimleri birbirinden bağımsız, birebir karşılaştırma
YAPILAMAZ. Karşılaştırma NİSPİ olmalı: "büyüme oranı referansın yüzde kaçına
düşüyor" ve "hangi şiddette tamamen infeasible oluyor" sorularının cevapları
karşılaştırılabilir, mutlak O2/glikoz sayıları değil.
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from extremofil_fba import bakim_reaksiyonunu_bul, mars_kisitlarini_uygula, modeli_yukle, referans_buyume

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SONUC_KLASORU = os.path.join(PROJE_KOKU, "results")
SOLVER_TOLERANCE = 1e-9

# Şiddet ekseni ankorları: t=0 çok sert (ama B. subtilis'in ilk varsayımı
# kadar keyfi değil -- ön taramada bu seviyede bile bakım x20'ye kadar
# feasible kaldığı görüldü), t=1 = modelin kendi kalibre edilmiş referansı.
SERT = dict(o2=0.02, glc=0.005, org=0.0005, h2o=0.05)
ILIMLI = dict(o2=3.15, glc=1.0, org=0.1, h2o=1000.0)

BAKIM_CARPANLARI = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
T_DEGERLERI = [round(i * 0.02, 2) for i in range(51)]


def kisit_degerleri(t):
    o2 = SERT["o2"] + t * (ILIMLI["o2"] - SERT["o2"])
    glc = SERT["glc"] + t * (ILIMLI["glc"] - SERT["glc"])
    org = SERT["org"] + t * (ILIMLI["org"] - SERT["org"])
    h2o = SERT["h2o"] + t * (ILIMLI["h2o"] - SERT["h2o"])
    return o2, glc, org, h2o


def tek_nokta_calistir(model, atpm, atpm_taban, t, bakim_carpani):
    o2, glc, org, h2o = kisit_degerleri(t)
    mars_kisitlarini_uygula(
        model, atpm, o2_ub=o2, glc_ub=glc, organik_ub=org, h2o_ub=h2o,
        bakim_carpani=bakim_carpani, bakim_taban=atpm_taban, sessiz=True,
    )
    sol = model.optimize()
    buyume = sol.objective_value if sol.status == "optimal" else None
    # solver toleransinin cok altindaki "teknik olarak optimal ama aslinda
    # sifir" sonuclari da sifir/infeasible olarak isaretle (bkz. docstring)
    if buyume is not None and buyume < SOLVER_TOLERANCE * 10:
        buyume = None
    return sol.status, buyume, o2, glc, org, h2o


def tarama_yap(model, atpm, atpm_taban):
    satirlar = []
    for bakim_x in BAKIM_CARPANLARI:
        for t in T_DEGERLERI:
            durum, buyume, o2, glc, org, h2o = tek_nokta_calistir(model, atpm, atpm_taban, t, bakim_x)
            satirlar.append(dict(
                bakim_carpani=bakim_x, t=t, O2_ub=o2, glc_ub=glc, organik_ub=org, h2o_ub=h2o,
                durum=durum, buyume=buyume,
            ))
    return pd.DataFrame(satirlar)


def grafik_ciz(df, baseline_buyume, dosya_yolu):
    fig, ax = plt.subplots(figsize=(8, 5))
    for bakim_x, grup in df.groupby("bakim_carpani"):
        gecerli = grup.dropna(subset=["buyume"]).sort_values("t")
        ax.plot(gecerli["t"], 100 * gecerli["buyume"] / baseline_buyume, marker=".", label=f"bakım ×{bakim_x}")
    ax.set_xlabel("Şiddet ekseni t  (0 = çok sert, 1 = model referansı)")
    ax.set_ylabel("Büyüme oranı (kendi referansına göre %)")
    ax.set_title("Salinibacter ruber (iMB631): Mars kısıt şiddetine duyarlılık")
    ax.legend(title="Bakım çarpanı")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(dosya_yolu, dpi=150)
    plt.close(fig)
    print(f"Grafik kaydedildi: {dosya_yolu}")


def main():
    os.makedirs(SONUC_KLASORU, exist_ok=True)

    model = modeli_yukle()
    model.solver.configuration.tolerances.feasibility = SOLVER_TOLERANCE
    baseline = referans_buyume(model)
    atpm = bakim_reaksiyonunu_bul(model)
    atpm_taban = atpm.upper_bound

    df = tarama_yap(model, atpm, atpm_taban)

    csv_yolu = os.path.join(SONUC_KLASORU, "extremofil_duyarlilik_sonuclari.csv")
    df.to_csv(csv_yolu, index=False)
    print(f"\nSonuçlar kaydedildi: {csv_yolu} ({len(df)} satır)")

    grafik_yolu = os.path.join(SONUC_KLASORU, "extremofil_buyume_vs_siddet.png")
    grafik_ciz(df, baseline.objective_value, grafik_yolu)

    print("\n--- Özet: her bakım çarpanı için ilk feasible nokta ---")
    for bakim_x, grup in df.groupby("bakim_carpani"):
        gecerli = grup.dropna(subset=["buyume"]).sort_values("t")
        if gecerli.empty:
            print(f"bakım ×{bakim_x}: taranan aralıkta hiçbir t değeri feasible değil")
            continue
        ilk = gecerli.iloc[0]
        pct = 100 * ilk["buyume"] / baseline.objective_value
        print(
            f"bakım ×{bakim_x}: ilk feasible t={ilk['t']:.2f} "
            f"-> büyüme={ilk['buyume']:.6f} (%{pct:.2f} referans)"
        )


if __name__ == "__main__":
    main()
