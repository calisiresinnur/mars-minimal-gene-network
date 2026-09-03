"""
GERÇEK minimal gen ağı inşası — Salinibacter ruber (iMB631), ardışık indirgeme.

Bkz. minimal_ag_insa.py (B. subtilis versiyonu) — aynı yöntem, aynı
gerekçe. Kullanıcının kritik sorusu üzerine: tekli gen silme testinin
"esansiyel gen" listesi, izoenzim/yedek-yol gruplarını (aynı reaksiyonu
yapabilen birden fazla gen) YAKALAYAMIYOR -- her biri tek başına
esansiyel değil ama HEPSİ birden çıkarılırsa ağ çöküyor. Bu script,
genleri TEK TEK, o ana kadar indirgenmiş ağa göre test ederek gerçekten
işlevsel bir minimal ağ inşa ediyor.
"""

import os

import pandas as pd

from extremofil_fba import bakim_reaksiyonunu_bul, mars_kisitlarini_uygula, modeli_yukle, referans_ortami_uygula
from extremofil_gen_silme import MARS_SENARYOLARI

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SONUC_KLASORU = os.path.join(PROJE_KOKU, "results")

ESANSIYELLIK_ESIGI = 0.01
SOLVER_TOLERANCE = 1e-9


def ardisik_indirgeme(model, esik=ESANSIYELLIK_ESIGI, sessiz=False):
    model.solver.configuration.tolerances.feasibility = SOLVER_TOLERANCE
    wt = model.optimize(raise_error=False)
    if model.solver.status != "optimal":
        raise RuntimeError(f"WT durumu optimal değil ({model.solver.status})")
    wt_buyume = wt.objective_value
    if not sessiz:
        print(f"WT büyüme: {wt_buyume:.6f}")

    tum_gen_id = [g.id for g in model.genes]
    cikarilan, tutulan = [], []

    for i, gid in enumerate(tum_gen_id):
        gene = model.genes.get_by_id(gid)
        if not gene.functional:
            continue
        # bkz. minimal_ag_insa.py -- "with model:" ile GEÇİCİ test, kalıcı
        # silme sadece bağlam DIŞINDA (gene.functional=True geri yükleme
        # yapmıyor, canlı yakalanan bir hata -- orada belgelendi).
        with model:
            gene.knock_out()
            sol = model.optimize(raise_error=False)
            buyume = sol.objective_value if model.solver.status == "optimal" else 0.0
        oran = (buyume / wt_buyume) if wt_buyume else 0.0
        if oran >= esik:
            gene.knock_out()
            cikarilan.append(gid)
        else:
            tutulan.append(gid)
        if not sessiz and (i + 1) % 100 == 0:
            print(f"  {i + 1}/{len(tum_gen_id)} -- tutulan: {len(tutulan)}, çıkarılan: {len(cikarilan)}")

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
            model, atpm, o2_ub=senaryo["o2"], glc_ub=senaryo["glc"], organik_ub=senaryo["org"],
            h2o_ub=senaryo["h2o"], bakim_carpani=senaryo["bakim_carpani"], bakim_taban=3.15, sessiz=True,
        )
    return uygula


def main():
    os.makedirs(SONUC_KLASORU, exist_ok=True)

    sonuclar = [senaryo_calistir("Referans", referans_ortami_uygula)]
    for s in MARS_SENARYOLARI:
        sonuclar.append(senaryo_calistir(s["etiket"], mars_kisiti(s)))

    eski = pd.read_csv(os.path.join(SONUC_KLASORU, "extremofil_gen_silme_sonuclari.csv"))

    print("\n\n=== ÖZET: her senaryo için tekli-silme vs gerçek minimal ağ ===")
    satirlar = []
    for etiket, tutulan, cikarilan, buyume in sonuclar:
        eski_esansiyel = set(eski[(eski.senaryo == etiket) & (eski.esansiyel == True)].gen_id)
        print(f"{etiket:26s}: tekli-silme={len(eski_esansiyel):4d}  gerçek_minimal={len(tutulan):4d}  "
              f"fark=+{len(tutulan - eski_esansiyel):3d}  büyüme={buyume}")
        satirlar.append(dict(senaryo=etiket, tekli_silme_sayisi=len(eski_esansiyel),
                              gercek_minimal_sayisi=len(tutulan), fark=len(tutulan - eski_esansiyel),
                              indirgenmis_ag_buyume=buyume))
        pd.DataFrame({"gen_id": sorted(tutulan), "durum": "tutuldu_gerekli"}).to_csv(
            os.path.join(SONUC_KLASORU, f"extremofil_minimal_ag_tutulan_genler_{etiket}.csv"), index=False)

    pd.DataFrame(satirlar).to_csv(os.path.join(SONUC_KLASORU, "extremofil_minimal_ag_ozet_tum_senaryolar.csv"), index=False)

    print("\n=== Referans minimal ağı ile Mars minimal ağları arasındaki fark ===")
    ref_tutulan = sonuclar[0][1]
    for etiket, tutulan, _, _ in sonuclar[1:]:
        sadece_mars = tutulan - ref_tutulan
        sadece_ref = ref_tutulan - tutulan
        print(f"{etiket}: Sadece Mars'ta gerekli: {len(sadece_mars)} gen "
              f"({', '.join(sorted(sadece_mars)) if sadece_mars else '-'}) | "
              f"Sadece referansta gerekli: {len(sadece_ref)} gen")

    print("\nKaydedildi: results/extremofil_minimal_ag_*.csv")


if __name__ == "__main__":
    main()
