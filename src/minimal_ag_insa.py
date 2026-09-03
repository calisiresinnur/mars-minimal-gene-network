"""
GERÇEK minimal gen ağı inşası — ardışık/açgözlü (greedy sequential) indirgeme.

## Neden bu script var (kullanıcının kritik sorusu üzerine, 2026-09-01)

`mars_gen_silme.py`'deki "esansiyel gen" listesi TEK TEK gen silme
testinden geliyor -- her gen, DİĞER HERKES modelde dururken tek başına
siliniyor. Bu genler "esansiyel" (WT durumda tek başına silinince
büyüme çöküyor) diye etiketleniyor, geri kalanlar "esansiyel değil".

AMA: "esansiyel olmayan TÜM genleri AYNI ANDA çıkarırsak ağ hâlâ çalışır
mı?" sorusu HİÇ TEST EDİLMEMİŞTİ. Test edilince (kullanıcının ısrarıyla)
sonuç ÇARPICI: 171 "esansiyel" genin oluşturduğu indirgenmiş ağ
İNFEASIBLE çıkıyor. Neden: glikoliz ve TCA döngüsünün TAMAMI (GLCpts,
PGI, PFK, FBA, GAPD, PGK, ENO, PYK, CS, ACONT, ICDHyr, AKGDH, SUCOAS,
FUM, MDH, PDH...) hiçbiri tek başına esansiyel değil (izoenzim/yedek
yolları var), ama HEPSİ AYNI ANDA çıkarılınca hücrenin enerji/karbon
üretecek HİÇBİR yolu kalmıyor. Tekli gen silme, bu tür "birlikte zorunlu
ama tek tek değil" (sinerjik/kombinatoryal esansiyellik) durumları
YAPISAL OLARAK yakalayamıyor.

## Doğru yöntem: ardışık/açgözlü indirgeme

`mars_gen_silme.py`'nin tek-seferlik "WT'ye göre esansiyel/değil" testi
yerine, genleri TEK TEK, HER SEFERİNDE O ANA KADAR İNDİRGENMİŞ AĞA göre
test ediyoruz:
  1. WT (tam model) büyümesini ölç.
  2. Genleri sırayla dolaş. Her gen için: O ANKİ (kısmen indirgenmiş)
     modelde bu geni sil, büyüme hâlâ eşiğin üzerinde mi kontrol et.
  3. Eşiğin üzerindeyse KALICI olarak sil (ağ küçülmeye devam eder).
     Değilse geri al (bu gen ağda kalır).
  4. Sona kadar devam et.

Bu yöntem izoenzim/yedek-yol durumlarını DOĞRU ele alıyor: bir izoenzim
çifti varsa, ilki silinince (o an hâlâ diğer izoenzim duruyor) sorun
çıkmaz, silinir. İkincisi silinmeye çalışılınca (o an İLK İZOENZİM ZATEN
YOK) büyüme çöker, bu yüzden ikincisi ağda TUTULUR. Sonuç: izoenzim
grubundan en az bir gen her zaman kalır -- tam da olması gereken.

NOT (dürüstlük): bu açgözlü/sıralı yöntem, genellerin test edilme
SIRASINA bağlı olarak farklı (ama hepsi geçerli/işlevsel) minimal
ağlara varabilir -- TEK bir "doğru" minimal ağ garantisi yok, ama
İŞLEVSEL bir minimal ağ garantisi var (ki bu, mars_gen_silme.py'nin
sağladığı "esansiyel gen listesi"nden -- ki işlevsel bile değil -- çok
daha güçlü bir iddia).
"""

import os
import sys

import pandas as pd

from mars_fba import bakim_reaksiyonunu_bul, mars_kisitlarini_uygula, modeli_yukle
from mars_gen_silme import MARS_SENARYOLARI

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SONUC_KLASORU = os.path.join(PROJE_KOKU, "results")

ESANSIYELLIK_ESIGI = 0.01
SOLVER_TOLERANCE = 1e-9


def ardisik_indirgeme(model, esik=ESANSIYELLIK_ESIGI, sessiz=False):
    """Modeli yerinde (in-place) indirger. Döner: (tutulan_genler, cikarilan_genler, son_buyume)."""
    model.solver.configuration.tolerances.feasibility = SOLVER_TOLERANCE
    wt = model.optimize(raise_error=False)
    if model.solver.status != "optimal":
        raise RuntimeError(f"WT durumu optimal değil ({model.solver.status}) -- indirgeme yapılamaz")
    wt_buyume = wt.objective_value
    if not sessiz:
        print(f"WT büyüme: {wt_buyume:.6f}")

    tum_gen_id = [g.id for g in model.genes]
    cikarilan = []
    tutulan = []

    for i, gid in enumerate(tum_gen_id):
        gene = model.genes.get_by_id(gid)
        if not gene.functional:
            # zaten (baska bir genin AND-baglantili GPR'i uzerinden) devre disi
            continue

        # DİKKAT (bu script'te canlı yakalanan bir hata): gene.functional = True
        # reaksiyon sınırlarını GERİ YÜKLEMİYOR (cobra'da bilinmeyen/dolaylı bir
        # davranış). Bunun yerine "with model:" bağlamı kullanılıyor -- bu,
        # GEÇİCİ test için bağlam çıkışında sınırları doğru şekilde eski haline
        # getiriyor (ayrıca doğrulandı). Kalıcı silme SADECE bağlam DIŞINDA yapılır.
        with model:
            gene.knock_out()
            sol = model.optimize(raise_error=False)
            buyume = sol.objective_value if model.solver.status == "optimal" else 0.0
        oran = (buyume / wt_buyume) if wt_buyume else 0.0
        if oran >= esik:
            gene.knock_out()  # KALICI -- bağlam dışında, ağ küçülmeye devam eder
            cikarilan.append(gid)
        else:
            tutulan.append(gid)  # dokunulmadı, ağda kalıyor
        if not sessiz and (i + 1) % 100 == 0:
            print(f"  {i + 1}/{len(tum_gen_id)} gen tarandi -- su ana kadar tutulan: {len(tutulan)}, cikarilan: {len(cikarilan)}")

    son_sol = model.optimize(raise_error=False)
    return tutulan, cikarilan, son_sol


def senaryo_calistir(etiket, kisit_uygula):
    model = modeli_yukle()
    if kisit_uygula is not None:
        kisit_uygula(model)
    print(f"\n=== {etiket} ===")
    tutulan, cikarilan, son_sol = ardisik_indirgeme(model)
    print(f"Toplam gen: {len(tutulan) + len(cikarilan)} | Tutulan: {len(tutulan)} | "
          f"Çıkarılan: {len(cikarilan)} | Son durum: {model.solver.status} | "
          f"Büyüme: {son_sol.objective_value}")
    return etiket, set(tutulan), set(cikarilan), son_sol.objective_value if model.solver.status == "optimal" else None


def mars_kisiti(senaryo):
    def uygula(model):
        atpm = bakim_reaksiyonunu_bul(model)
        mars_kisitlarini_uygula(
            model, atpm, o2_lb=senaryo["o2"], glc_lb=senaryo["glc"], h2o_cap=senaryo["h2o"],
            bakim_carpani=senaryo["bakim_carpani"], sessiz=True,
        )
    return uygula


def main():
    os.makedirs(SONUC_KLASORU, exist_ok=True)

    sonuclar = [senaryo_calistir("Dunya_referans", None)]
    for s in MARS_SENARYOLARI:
        sonuclar.append(senaryo_calistir(s["etiket"], mars_kisiti(s)))

    # Karsilastirma: mars_gen_silme.py'nin (yanlis) "esansiyel gen" sayisiyla
    eski = pd.read_csv(os.path.join(SONUC_KLASORU, "gen_silme_sonuclari.csv"))

    print("\n\n=== ÖZET: her senaryo için tekli-silme vs gerçek minimal ağ ===")
    satirlar = []
    for etiket, tutulan, cikarilan, buyume in sonuclar:
        eski_esansiyel = set(eski[(eski.senaryo == etiket) & (eski.esansiyel == True)].gen_id)
        print(f"{etiket:32s}: tekli-silme={len(eski_esansiyel):4d}  gerçek_minimal={len(tutulan):4d}  "
              f"fark=+{len(tutulan - eski_esansiyel):3d}  büyüme={buyume}")
        satirlar.append(dict(senaryo=etiket, tekli_silme_sayisi=len(eski_esansiyel),
                              gercek_minimal_sayisi=len(tutulan), fark=len(tutulan - eski_esansiyel),
                              indirgenmis_ag_buyume=buyume))
        pd.DataFrame({"gen_id": sorted(tutulan), "durum": "tutuldu_gerekli"}).to_csv(
            os.path.join(SONUC_KLASORU, f"minimal_ag_tutulan_genler_{etiket}.csv"), index=False)

    pd.DataFrame(satirlar).to_csv(os.path.join(SONUC_KLASORU, "minimal_ag_ozet_tum_senaryolar.csv"), index=False)

    # Dunya minimal ag ile her Mars minimal agini karsilastir
    print("\n=== Dünya minimal ağı ile Mars minimal ağları arasındaki fark ===")
    dunya_tutulan = sonuclar[0][1]
    for etiket, tutulan, _, _ in sonuclar[1:]:
        sadece_mars = tutulan - dunya_tutulan
        sadece_dunya = dunya_tutulan - tutulan
        print(f"{etiket}: Sadece Mars'ta gerekli: {len(sadece_mars)} gen "
              f"({', '.join(sorted(sadece_mars)) if sadece_mars else '-'}) | "
              f"Sadece Dünya'da gerekli: {len(sadece_dunya)} gen")

    print("\nKaydedildi: results/minimal_ag_tutulan_genler_<senaryo>.csv, "
          "results/minimal_ag_ozet_tum_senaryolar.csv")


if __name__ == "__main__":
    main()
