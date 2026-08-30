"""
Karşılaştırmalı genomik: FBA gen esansiyellik sonuçlarını gerçek B. subtilis
stres-regulonu / DEG (differentially expressed genes) literatürüyle karşılaştırma.

VERİ KAYNAĞI (literatür taraması, 2026-08-30):
data/stres_regulon_genleri.csv dosyasındaki gen listesi 4 gerçek kaynaktan
derlendi:
  - PerR regulonu (oksidatif stres): Fuangthong ve ark. 2002, J Bacteriol
    (PMID 12029044). PROTECT/EXPOSE-E deneyinde (Moeller ve ark. 2012,
    B. subtilis sporları 559 gün gerçek uzay + simüle Mars koşullarına
    maruz bırakıldı) PerR regulonu Mars-simule koşullarda indüklenmiş bulundu.
  - SigV regulonu (hücre zarfı/lizozim direnci): Guariglia-Oropeza & Helmann
    2011, J Bacteriol (PMID 21926231).
  - ResD-ResE regulonu (O2 kısıtlaması, anaerobik solunum/fermantasyon):
    Nakano laboratuvarının klasik çalışmaları (resABCDE, cydABCD, nasBCDEF,
    fnr vb.).
  - ISS uçuş deneyi (BRIC-21/BRIC-23): Nickerson lab, npj Microgravity 2019
    (PMC6323116) -- iki ayrı ISS misyonunda TUTARLI bulunan 91 DEG'den O2
    kısıtlamasıyla ilişkili olanlar (yer kontrolünde/O2-zengin ortamda
    indüklenen anaerobik/fermantasyon genleri).

ÖNEMLİ KISIT: iYO844 SADECE metabolik (enzim kodlayan) genleri içeriyor --
sigma faktörleri, transkripsiyon faktörleri ve düzenleyici proteinler (ör.
PerR, Fur, ResD, ResE, LexA, RecA, CtsR, SigV, Clp proteaz sistemi) modelde
YOK çünkü bunlar doğrudan bir metabolik reaksiyon katalizlemiyor. Bu, GEM'in
kapsamının doğal bir sınırı -- 68 aday gen isminden sadece 29'u modelde
bulundu (bkz. data/stres_regulon_genleri.csv'nin nasıl derlendiği).

BULGU: Aşağıdaki cross-reference, mars_gen_silme.py'nin bulduğu 3 genin
(pabB, menC, menD) rastgele değil, gerçek O2-kısıtlaması biyolojisiyle
tutarlı bir yerde durduğunu gösteriyor: menC/menD PRİMER solunum zinciri
kofaktörü (menakinon) biyosentezinde -- Mars'ta O2 zaten o kadar kısıtlı ki
bu yol darboğaz değil. Buna karşılık, gerçek B. subtilis biyolojisinde
BİLİNEN alternatif/yedek O2-kısıtlaması yolları (cydAB, nasBCDEF, narGHJK,
fermantasyon genleri ldh/lctP/bdhA) modelin HİÇBİR senaryosunda esansiyel
çıkmıyor -- beklenen ve tutarlı bir sonuç, çünkü "alternatif/yedek yol"
olmaları zaten esansiyel olmamaları gerektiği anlamına geliyor.

Çıktı: results/deg_karsilastirma.csv
"""

import os

import pandas as pd

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SONUC_KLASORU = os.path.join(PROJE_KOKU, "results")
STRES_GEN_CSV = os.path.join(PROJE_KOKU, "data", "stres_regulon_genleri.csv")
GEN_SILME_CSV = os.path.join(SONUC_KLASORU, "gen_silme_sonuclari.csv")


def main():
    if not os.path.exists(GEN_SILME_CSV):
        raise FileNotFoundError(
            f"{GEN_SILME_CSV} bulunamadı -- önce `python src/mars_gen_silme.py` çalıştırılmalı."
        )

    stres_genleri = pd.read_csv(STRES_GEN_CSV)
    gen_silme = pd.read_csv(GEN_SILME_CSV)
    pivot = gen_silme.pivot(index="gen_id", columns="senaryo", values="esansiyel")

    birlesik = stres_genleri.merge(pivot, left_on="gen_id", right_index=True, how="left")
    senaryo_kolonlari = [c for c in pivot.columns]
    modelde_yok = birlesik[birlesik[senaryo_kolonlari[0]].isna()]
    if not modelde_yok.empty:
        print(f"UYARI: {len(modelde_yok)} gen data/stres_regulon_genleri.csv içinde ama "
              f"gen_silme_sonuclari.csv'de yok -- önce mars_gen_silme.py çalıştırıldığından emin olun.")

    csv_yolu = os.path.join(SONUC_KLASORU, "deg_karsilastirma.csv")
    birlesik.to_csv(csv_yolu, index=False)

    print(f"Karşılaştırma tablosu kaydedildi: {csv_yolu} ({len(birlesik)} gen)\n")
    pd.set_option("display.width", 200)
    print(birlesik[["gen_id", "isim", "kategori"] + senaryo_kolonlari].to_string(index=False))

    print("\n--- Özet: kategoriye göre kaç gen esansiyel ---")
    for kategori, grup in birlesik.groupby("kategori"):
        sayilar = ", ".join(f"{k}={int(grup[k].sum())}" for k in senaryo_kolonlari)
        print(f"{kategori:55s} ({len(grup)} gen): {sayilar}")


if __name__ == "__main__":
    main()
