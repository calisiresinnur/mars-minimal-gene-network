"""
Mars kısıt şiddetine duyarlılık analizi (sensitivity analysis).

NEDEN BU SCRIPT VAR (literatür taramasından sonra eklendi):
"Kanonik" tek bir Mars senaryosu (belirli bir O2/glikoz/su/bakım kombinasyonu)
iddia etmek yerine bilinçli olarak bundan vazgeçtik. Sebep:

  1. Radyasyonun ATP bakım maliyetine (ATPM) etkisini sayısallaştıran bir
     kaynak yok. Genome-scale modellerin astrobiyoloji uygulamalarını
     inceleyen güncel bir derleme (Noirungsee ve ark. 2024, Environ Microbiol
     Reports, PMC10866088) bile "radyasyon etkileri analize dahil edilebilir"
     diyor ama somut bir sayısal yöntem önermiyor.
  2. Gerçek Mars atmosfer yüzdelerini (Cortesão ve ark. 2019, Front Microbiol,
     PMC6399134: %95.54 CO2, %0.13 O2, %0.03 H2O,
     toplam basın 0.69 kPa) doğrudan bir FBA akı sınırına (mmol/gDW/h)
     çevirecek kinetik/taşınım (mass-transfer) verisi mevcut değil -- FBA bir
     sabit-akı-sınırı formalizmi, kısmi basınç oranını orantılı bir akı
     sınırına çevirmek bilimsel olarak savunulamaz bir kesinlik iddiası olur.

Bu yüzden doğru ve savunulabilir yaklaşım: kısıt şiddetini TEK BİR ortak
eksende (t: 0=en sert ilk varsayım, 1=daha ılımlı bir uç) tarayıp, farklı
bakım-enerjisi çarpanları (radyasyon onarım maliyeti varsayımı, kendisi de
bilinmeyen bir parametre) için büyüme oranının bu şiddete nasıl tepki
verdiğini raporlamak. Makalenin Sonuçlar bölümünde "tek sayı" yerine bu
duyarlılık eğrileri sunulmalı -- bu, mars_kalibrasyon.py'daki nokta
taramasının sistematik/görselleştirilmiş hali.

Kaynaklar için bkz. README.md > Kaynaklar.
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from mars_fba import bakim_reaksiyonunu_bul, mars_kisitlarini_uygula, modeli_yukle, referans_buyume

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SONUC_KLASORU = os.path.join(PROJE_KOKU, "results")

# Şiddet ekseni ankorları: t=0 en sert ilk varsayım (mars_fba.py'nin öntanımlısı),
# t=1 daha ılımlı bir uç. Aradaki her t için üç kısıt BİRLİKTE doğrusal olarak
# gevşetiliyor -- bunlar bağımsız/literatürden gelen kesin sayılar değil,
# şiddet ekseninin uç noktalarını tanımlayan keyfi ama makul ankorlar.
SERT = dict(o2=-0.5, glc=-0.05, h2o=1.0)
ILIMLI = dict(o2=-10.0, glc=-2.0, h2o=30.0)

BAKIM_CARPANLARI = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
T_DEGERLERI = [round(i * 0.02, 2) for i in range(51)]  # 0.00, 0.02, ..., 1.00


def kisit_degerleri(t):
    o2 = SERT["o2"] + t * (ILIMLI["o2"] - SERT["o2"])
    glc = SERT["glc"] + t * (ILIMLI["glc"] - SERT["glc"])
    h2o = SERT["h2o"] + t * (ILIMLI["h2o"] - SERT["h2o"])
    return o2, glc, h2o


def tek_nokta_calistir(model, atpm, atpm_taban, t, bakim_carpani):
    o2, glc, h2o = kisit_degerleri(t)
    # bakim_taban'ı MUTLAKA açıkça veriyoruz: model/atpm bu döngü boyunca tekrar
    # kullanılıyor, atpm.lower_bound bir önceki çağrıdan kalma mutasyona uğramış
    # bir değer -- taban olarak kullanılırsa katlanarak büyür (bkz. mars_fba.py
    # içindeki uyarı). atpm_taban, döngü başında modelin ORİJİNAL değerinden
    # bir kere yakalanmış sabit bir sayı.
    mars_kisitlarini_uygula(
        model, atpm, o2_lb=o2, glc_lb=glc, h2o_cap=h2o, bakim_carpani=bakim_carpani,
        bakim_taban=atpm_taban, sessiz=True,
    )
    sol = model.optimize()
    return model.solver.status, sol.objective_value, o2, glc, h2o


def tarama_yap(model, atpm, atpm_taban):
    satirlar = []
    for bakim_x in BAKIM_CARPANLARI:
        for t in T_DEGERLERI:
            durum, buyume, o2, glc, h2o = tek_nokta_calistir(model, atpm, atpm_taban, t, bakim_x)
            satirlar.append(
                dict(
                    bakim_carpani=bakim_x,
                    t=t,
                    O2_lb=o2,
                    glc_lb=glc,
                    h2o_cap=h2o,
                    durum=durum,
                    buyume=buyume if durum == "optimal" else None,
                )
            )
    return pd.DataFrame(satirlar)


def grafik_ciz(df, baseline_buyume, dosya_yolu):
    fig, ax = plt.subplots(figsize=(8, 5))
    for bakim_x, grup in df.groupby("bakim_carpani"):
        gecerli = grup.dropna(subset=["buyume"]).sort_values("t")
        ax.plot(gecerli["t"], 100 * gecerli["buyume"] / baseline_buyume, marker=".", label=f"bakım ×{bakim_x}")
    ax.set_xlabel("Şiddet ekseni t  (0 = en sert ilk varsayım, 1 = ılımlı uç)")
    ax.set_ylabel("Büyüme oranı (Dünya benzeri referansa göre %)")
    ax.set_title("Mars kısıt şiddetine duyarlılık analizi")
    ax.legend(title="Bakım çarpanı")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(dosya_yolu, dpi=150)
    plt.close(fig)
    print(f"Grafik kaydedildi: {dosya_yolu}")


def main():
    os.makedirs(SONUC_KLASORU, exist_ok=True)

    model = modeli_yukle()
    baseline = referans_buyume(model)  # kısıtlar uygulanmadan ÖNCE ölçülmeli
    atpm = bakim_reaksiyonunu_bul(model)
    atpm_taban = atpm.lower_bound  # orijinal 9.0 -- döngü boyunca hep bundan çarpıyoruz

    df = tarama_yap(model, atpm, atpm_taban)

    csv_yolu = os.path.join(SONUC_KLASORU, "duyarlilik_sonuclari.csv")
    df.to_csv(csv_yolu, index=False)
    print(f"\nSonuçlar kaydedildi: {csv_yolu} ({len(df)} satır)")

    grafik_yolu = os.path.join(SONUC_KLASORU, "buyume_vs_siddet.png")
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
            f"(O2={ilk['O2_lb']:.2f}, glc={ilk['glc_lb']:.2f}, h2o=±{ilk['h2o_cap']:.2f}) "
            f"-> büyüme={ilk['buyume']:.4f} (%{pct:.1f} Dünya)"
        )


if __name__ == "__main__":
    main()
