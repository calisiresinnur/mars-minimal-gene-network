"""
Ardışık/açgözlü minimal ağ indirgemesinin SIRA-SAĞLAMLIĞI testi.

NEDEN BU SCRIPT VAR: `minimal_ag_insa.py`'nin ilk çalıştırmasında "Mars
minimal ağı Dünya'dan 14-16 gen daha büyük" denmişti — TEK bir gen
sıralamasına (modelin kendi varsayılan gen listeleme sırası) dayanarak.
Ardışık indirgeme yöntemi, izoenzim gruplarından HANGİSİNİN tutulacağı
konusunda test SIRASINA duyarlı olabilir (bkz. minimal_ag_insa.py'nin
docstring'indeki dürüstlük notu). Bu script, farklı sıralamalarla
(orijinal, tam ters, 4 rastgele) aynı karşılaştırmayı tekrarlayıp
bulgunun ne kadarının sağlam (sıradan bağımsız), ne kadarının sıraya
duyarlı olduğunu ölçüyor.

SONUÇ (2026-09-01): 6 denemede fark (Mars−Dünya) −4 ile +11 arasında,
5/6'sı pozitif — YÖN genel olarak sağlam ama BÜYÜKLÜK değil. Solunum
elektron taşıma zinciri (CYOO3/CYOR3m) iki farklı rastgele sırada da
BİREBİR aynı genlerle Mars'a özgü çıktı — bu tema sağlam. Detay:
README > "Mars senaryoları da yeniden test edildi".
"""

import random

import mars_fba as mf
from mars_gen_silme import MARS_SENARYOLARI

SOLVER_TOLERANCE = 1e-9
ESANSIYELLIK_ESIGI = 0.01


def ardisik_indirgeme_sira(model, sira, esik=ESANSIYELLIK_ESIGI):
    model.solver.configuration.tolerances.feasibility = SOLVER_TOLERANCE
    wt = model.optimize(raise_error=False)
    wt_buyume = wt.objective_value
    tutulan = []
    for gid in sira:
        gene = model.genes.get_by_id(gid)
        if not gene.functional:
            continue
        with model:
            gene.knock_out()
            sol = model.optimize(raise_error=False)
            buyume = sol.objective_value if model.solver.status == "optimal" else 0.0
        oran = (buyume / wt_buyume) if wt_buyume else 0.0
        if oran >= esik:
            gene.knock_out()
        else:
            tutulan.append(gid)
    return set(tutulan)


def dunya_mars_farki(sira, mars_senaryo=MARS_SENARYOLARI[0]):
    model_e = mf.modeli_yukle()
    tutulan_e = ardisik_indirgeme_sira(model_e, sira)

    model_m = mf.modeli_yukle()
    atpm = mf.bakim_reaksiyonunu_bul(model_m)
    mf.mars_kisitlarini_uygula(
        model_m, atpm, o2_lb=mars_senaryo["o2"], glc_lb=mars_senaryo["glc"],
        h2o_cap=mars_senaryo["h2o"], bakim_carpani=mars_senaryo["bakim_carpani"], sessiz=True,
    )
    tutulan_m = ardisik_indirgeme_sira(model_m, sira)

    return tutulan_e, tutulan_m


def main():
    model_ref = mf.modeli_yukle()
    tum_genler = [g.id for g in model_ref.genes]

    denemeler = [("orijinal", tum_genler)]
    denemeler.append(("ters", tum_genler[::-1]))
    random.seed(42)
    for i in range(4):
        sira = tum_genler[:]
        random.shuffle(sira)
        denemeler.append((f"rastgele_{i}", sira))

    print("=== Sıra-sağlamlık testi: B. subtilis, Dünya vs Mars_bakim_x1.5 ===\n")
    farklar = []
    for etiket, sira in denemeler:
        tutulan_e, tutulan_m = dunya_mars_farki(sira)
        fark = len(tutulan_m) - len(tutulan_e)
        sadece_mars = tutulan_m - tutulan_e
        print(f"{etiket:12s}: Dünya={len(tutulan_e):4d}  Mars={len(tutulan_m):4d}  "
              f"fark={fark:+3d}  sadece_Mars_genleri={sorted(sadece_mars)}")
        farklar.append(fark)

    print(f"\nTüm farklar: {farklar}")
    print(f"Pozitif (Mars≥Dünya) oranı: {sum(1 for f in farklar if f >= 0)}/{len(farklar)}")


if __name__ == "__main__":
    main()
