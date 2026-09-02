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

from extremofil_fba import modeli_yukle, referans_ortami_uygula

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


def main():
    os.makedirs(SONUC_KLASORU, exist_ok=True)

    print("=== Salinibacter ruber (iMB631) -- kalibre edilmiş referans için ARDIŞIK indirgeme ===")
    model = modeli_yukle()
    referans_ortami_uygula(model)
    tutulan, cikarilan, son_sol = ardisik_indirgeme(model)

    print(f"\nToplam gen: {len(tutulan) + len(cikarilan)}")
    print(f"Tutulan (gerçekten gerekli) gen: {len(tutulan)}")
    print(f"Çıkarılan (gerçekten gereksiz) gen: {len(cikarilan)}")
    print(f"İndirgenmiş ağ son durum: {model.solver.status}, büyüme: {son_sol.objective_value}")

    eski_yol = os.path.join(SONUC_KLASORU, "extremofil_gen_silme_sonuclari.csv")
    if os.path.exists(eski_yol):
        eski = pd.read_csv(eski_yol)
        kolon_adi = "senaryo" if "senaryo" in eski.columns else None
        if kolon_adi:
            ilk_senaryo = eski[kolon_adi].unique()[0]
            eski_esansiyel = set(eski[(eski[kolon_adi] == ilk_senaryo) & (eski.esansiyel == True)].gen_id)
            print(f"\nKarşılaştırma: tekli-silme esansiyel sayısı ({ilk_senaryo}): {len(eski_esansiyel)}")
            print(f"Bu script'in tuttuğu gen sayısı: {len(tutulan)}")

    pd.DataFrame({"gen_id": tutulan, "durum": "tutuldu_gerekli"}).to_csv(
        os.path.join(SONUC_KLASORU, "extremofil_minimal_ag_tutulan_genler.csv"), index=False)
    pd.DataFrame({"gen_id": cikarilan, "durum": "cikarildi_gereksiz"}).to_csv(
        os.path.join(SONUC_KLASORU, "extremofil_minimal_ag_cikarilan_genler.csv"), index=False)
    print("\nKaydedildi: results/extremofil_minimal_ag_*.csv")


if __name__ == "__main__":
    main()
