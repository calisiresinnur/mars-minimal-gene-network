"""
Tekli gen silme (single gene deletion) analizi.

Her senaryo için (Dünya benzeri referans + duyarlılık analizindeki 3 bakım
çarpanının "hafif marj" noktası) modelin 844 geninin her biri tek tek
"silinip" (ilgili reaksiyonlar kapatılıp) büyüme oranı yeniden hesaplanıyor.
Bir gen, o senaryoda WT (silinmemiş) büyümenin ESANSIYELLIK_ESIGI'nden azına
düşürüyorsa "esansiyel" sayılıyor.

NEDEN "sıkı marj" (t*+0.01) ve DAHA GEVŞEK NOKTALAR DEĞİL:
İlk denemede "hafif marj" (t*+0.05, büyüme ~%11.5 Dünya) kullanıldı ve
Dünya'yla TAMAMEN AYNI 171 esansiyel gen bulundu -- hiçbir Mars'a özgü etki
yoktu. Tam sınırda (t*) ise WT büyüme ~0'a çok yakın olduğu için esansiyellik
oranı (KO_büyüme / WT_büyüme) sayısal olarak anlamsızlaşır. Bu ikisi arasında,
sınıra ÇOK yakın bir nokta (t*+0.01, büyüme ~%2-23 Dünya, senaryoya göre
değişir) denendiğinde fark ortaya çıktı: 4 gen (folEA, pabB, menC, menD)
Mars'ta esansiyel OLMAKTAN ÇIKIYOR (tam tersi beklenen yönde -- yeni esansiyel
gen değil, esansiyelliğini kaybeden gen). Yorum: bu genler menakinon
(solunum zinciri kofaktörü) ve folat/tek-karbon metabolizması yollarına
hizmet ediyor; O2/karbon zaten o kadar kısıtlı ki darboğaz bu yollara hiç
ulaşmıyor, dolayısıyla bu genler "gerekli" değil "gereksiz" hale geliyor.
Bu etki SADECE sınıra çok yakın noktalarda görünüyor -- bkz. README >
Tekli gen silme bulguları.

Çıktı:
  - results/gen_silme_sonuclari.csv          (senaryo x gen x büyüme, ham veri)
  - results/mars_ozel_esansiyel_genler.csv   (SADECE Mars senaryolarında
    esansiyel, Dünya'da esansiyel OLMAYAN genler -- makalenin ana bulgu adayı)

KRİTİK DÜZELTME (kullanıcının bilimsel doğruluk denetimi sırasında bulundu,
2026-08-30): İlk çalıştırmada solver'ın VARSAYILAN feasibility tolerance'ı
(1e-7) kullanılmıştı. "Sıkı marj" senaryolarında WT büyüme zaten çok küçük
olduğundan (~0.0027/saat), biyokütle denklemindeki bazı kofaktörlerin
(ör. menakinon, mql7_c) gerektirdiği akı da çok küçük çıkıyor -- bazı
durumlarda solver toleransının sadece ~7 katı büyüklükte (7.2e-7 vs 1e-7).
Bu, gen silindiğinde LP çözücünün gerçekte imkânsız olan bir akıyı "toleransa
sığıyor" diye feasible/sıfırdan farklı büyüme olarak kabul etmesine yol açtı
-- yani "pabB/menC/menD Mars'ta esansiyellikten çıkıyor" bulgusu SAYISAL BİR
ARTEFAKTTI, gerçek bir biyolojik/model bulgusu değildi. Tolerance 1e-9'a
çekilince (aşağıdaki SOLVER_TOLERANCE), bu 4 genin (pabB, menC, menD, folEA)
DÖRDÜ DE Mars'ın üç senaryosunda da Dünya'daki gibi TAM ESANSİYEL çıkıyor --
yani düzeltilmiş bulgu: hiçbir gen esansiyellik durumunu değiştirmiyor.
"""

import os

import pandas as pd
from cobra.flux_analysis import single_gene_deletion

from mars_fba import bakim_reaksiyonunu_bul, mars_kisitlarini_uygula, modeli_yukle

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SONUC_KLASORU = os.path.join(PROJE_KOKU, "results")

ESANSIYELLIK_ESIGI = 0.01  # KO büyüme / WT büyüme bu değerin altındaysa "esansiyel"
ISLEMCI_SAYISI = 4  # bkz. README; makinede daha fazla çekirdek varsa artırılabilir
SOLVER_TOLERANCE = 1e-9  # varsayılan 1e-7 -- Mars senaryolarında WT büyüme çok küçük
# olduğu için gerekli kofaktör akışları toleransa yakın kalıyor ve gevşek tolerans
# sahte-feasible sonuçlar üretiyor (bkz. yukarıdaki KRİTİK DÜZELTME notu)

# mars_duyarlilik.py'deki bisection ile bulunan t* (tam feasibility sınırı)
# değerlerinin 0.01 ötesi -- yani "sınıra çok yakın, ama numerik olarak
# kullanılabilir bir WT büyümesi olan" noktalar (bkz. yukarıdaki NEDEN notu).
MARS_SENARYOLARI = [
    dict(etiket="Mars_bakim_x1.5_siki_marj", o2=-3.3861, glc=-0.6424, h2o=9.8102, bakim_carpani=1.5),
    dict(etiket="Mars_bakim_x2.0_siki_marj", o2=-4.4643, glc=-0.8637, h2o=13.1017, bakim_carpani=2.0),
    dict(etiket="Mars_bakim_x3.0_siki_marj", o2=-6.6208, glc=-1.3064, h2o=19.6847, bakim_carpani=3.0),
]


def senaryo_calistir(etiket, kisit_uygula):
    """kisit_uygula(model) -> None; None geçilirse Dünya benzeri referans (kısıtsız)."""
    model = modeli_yukle()
    if kisit_uygula is not None:
        kisit_uygula(model)
    model.solver.configuration.tolerances.feasibility = SOLVER_TOLERANCE

    wt = model.optimize()
    print(f"{etiket}: WT büyüme = {wt.objective_value:.6f} (durum: {model.solver.status})")

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
            model, atpm,
            o2_lb=senaryo["o2"], glc_lb=senaryo["glc"], h2o_cap=senaryo["h2o"],
            bakim_carpani=senaryo["bakim_carpani"], sessiz=True,
        )
    return uygula


def main():
    os.makedirs(SONUC_KLASORU, exist_ok=True)

    tum_sonuclar = [senaryo_calistir("Dunya_referans", None)]
    for s in MARS_SENARYOLARI:
        tum_sonuclar.append(senaryo_calistir(s["etiket"], mars_kisiti(s)))

    df = pd.concat(tum_sonuclar, ignore_index=True)
    csv_yolu = os.path.join(SONUC_KLASORU, "gen_silme_sonuclari.csv")
    df.to_csv(csv_yolu, index=False)
    print(f"\nHam sonuçlar kaydedildi: {csv_yolu} ({len(df)} satır)")

    pivot = df.pivot(index="gen_id", columns="senaryo", values="esansiyel")
    dunya_esansiyel = pivot["Dunya_referans"]
    mars_kolonlari = [c for c in pivot.columns if c != "Dunya_referans"]

    # Yön 1: Dünya'da esansiyel OLMAYAN ama en az bir Mars senaryosunda esansiyel olan genler
    # (beklenen yön: Mars YENİ esansiyel genler yaratır)
    mars_yeni_esansiyel = pivot[(~dunya_esansiyel) & pivot[mars_kolonlari].any(axis=1)].copy()
    mars_yeni_esansiyel["kac_mars_senaryosunda"] = mars_yeni_esansiyel[mars_kolonlari].sum(axis=1)
    mars_yeni_esansiyel = mars_yeni_esansiyel.sort_values("kac_mars_senaryosunda", ascending=False)
    mars_yeni_esansiyel.to_csv(os.path.join(SONUC_KLASORU, "mars_yeni_esansiyel_genler.csv"))

    # Yön 2: Dünya'da esansiyel OLAN ama en az bir Mars senaryosunda esansiyel OLMAYAN genler
    # (gözlenen asıl bulgu -- bkz. README > Tekli gen silme bulguları)
    mars_dispanse = pivot[dunya_esansiyel & (~pivot[mars_kolonlari]).any(axis=1)].copy()
    mars_dispanse["kac_mars_senaryosunda_dispanse"] = (~mars_dispanse[mars_kolonlari]).sum(axis=1)
    mars_dispanse = mars_dispanse.sort_values("kac_mars_senaryosunda_dispanse", ascending=False)
    mars_dispanse.to_csv(os.path.join(SONUC_KLASORU, "mars_dispanse_olan_genler.csv"))

    print("\n--- Ozet ---")
    for kolon in pivot.columns:
        print(f"{kolon:32s}: {int(pivot[kolon].sum()):4d} esansiyel gen / {len(pivot)}")
    print(f"\nMars'a ozgu YENI esansiyel gen adayi (beklenen yon): {len(mars_yeni_esansiyel)}")
    print(f"Mars'ta esansiyel OLMAKTAN CIKAN gen (gozlenen bulgu): {len(mars_dispanse)}")
    if len(mars_dispanse) > 0:
        print("  Bu genler:", ", ".join(mars_dispanse.index.tolist()))
    print(f"Kaydedildi: results/mars_yeni_esansiyel_genler.csv, results/mars_dispanse_olan_genler.csv")


if __name__ == "__main__":
    main()
