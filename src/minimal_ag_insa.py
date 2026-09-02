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


def main():
    os.makedirs(SONUC_KLASORU, exist_ok=True)

    print("=== B. subtilis (iYO844) -- Dünya benzeri referans için ARDIŞIK indirgeme ===")
    model = modeli_yukle()
    tutulan, cikarilan, son_sol = ardisik_indirgeme(model)

    print(f"\nToplam gen: {len(tutulan) + len(cikarilan)}")
    print(f"Tutulan (gerçekten gerekli) gen: {len(tutulan)}")
    print(f"Çıkarılan (gerçekten gereksiz) gen: {len(cikarilan)}")
    print(f"İndirgenmiş ağ son durum: {model.solver.status}, büyüme: {son_sol.objective_value}")

    # Karsilastirma: mars_gen_silme.py'nin (yanlis) "esansiyel gen" sayisiyla
    eski = pd.read_csv(os.path.join(SONUC_KLASORU, "gen_silme_sonuclari.csv"))
    eski_esansiyel = set(eski[(eski.senaryo == "Dunya_referans") & (eski.esansiyel == True)].gen_id)
    print(f"\nKarşılaştırma: mars_gen_silme.py'nin (tekli-silme) esansiyel sayısı: {len(eski_esansiyel)}")
    print(f"Bu script'in (ardışık, gerçekten işlevsel) tuttuğu gen sayısı: {len(tutulan)}")
    fark = set(tutulan) - eski_esansiyel
    print(f"Sadece bu yöntemde 'gerekli' çıkan ama tekli-silmede 'esansiyel değil' denen gen sayısı: {len(fark)}")

    pd.DataFrame({"gen_id": tutulan, "durum": "tutuldu_gerekli"}).to_csv(
        os.path.join(SONUC_KLASORU, "minimal_ag_tutulan_genler_Dunya.csv"), index=False)
    pd.DataFrame({"gen_id": cikarilan, "durum": "cikarildi_gereksiz"}).to_csv(
        os.path.join(SONUC_KLASORU, "minimal_ag_cikarilan_genler_Dunya.csv"), index=False)
    print("\nKaydedildi: results/minimal_ag_tutulan_genler_Dunya.csv, results/minimal_ag_cikarilan_genler_Dunya.csv")


if __name__ == "__main__":
    main()
