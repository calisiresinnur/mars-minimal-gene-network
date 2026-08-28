"""
Mars Yüzey Koşulları İçin Minimal Gen Ağı Modellemesi
IAC 2026 - Paper ID 114761

Bu script:
  1) Kürasyonu yapılmış bir genom-ölçekli metabolik modeli (iYO844, B. subtilis) yükler
  2) Referans (Dünya benzeri) büyüme oranını hesaplar
  3) Mars yüzey koşullarını (O2, CO2, organik karbon, su, radyasyon->bakım enerjisi)
     sayısal kısıtlara çevirip modele uygular
  4) Mars koşulunda büyüme oranını yeniden hesaplar ve iki sonucu karşılaştırır

Not: Buradaki kısıt değerleri (-0.5, -0.05, x3 vb.) kesin ölçümler değil, literatürden
türetilmiş gerekçeli ilk varsayımlardır. Makalenin Yöntem/Sınırlamalar bölümünde bu
şekilde açıkça ifade edilecek ve ileride bir duyarlılık analizi ile sağlamlaştırılabilir.
"""

import cobra


def modeli_yukle():
    model = cobra.io.load_model("iYO844")
    print(f"Model yüklendi: {len(model.reactions)} reaksiyon, {len(model.genes)} gen")
    return model


def referans_buyume(model):
    baseline = model.optimize()
    print(f"Dünya benzeri büyüme oranı (1/saat): {baseline.objective_value}")
    return baseline


def bakim_reaksiyonunu_bul(model):
    try:
        atpm = model.reactions.get_by_id("ATPM")
    except KeyError:
        atpm = [r for r in model.reactions if "maintenance" in r.name.lower()][0]
    print(f"Bakım reaksiyonu: {atpm.id} | mevcut sınırlar: {atpm.bounds}")
    return atpm


def mars_kisitlarini_uygula(model, atpm, bakim_carpani=3):
    # O2: Mars atmosferinin sadece ~%0.17'si O2 + toplam basınç Dünya'nın ~%0.6'sı -> ciddi kısıtla
    model.reactions.EX_o2_e.lower_bound = -0.5

    # CO2: Mars atmosferinin ~%95'i CO2 -> alıma açıkça izin ver, kısıtlayıcı olmasın
    model.reactions.EX_co2_e.bounds = (-1000, 1000)

    # Organik karbon: Mars yüzeyinde serbest glikoz/organik karbon yok denecek kadar az
    model.reactions.EX_glc__D_e.lower_bound = -0.05

    # Su: düşük su aktivitesi (aw ~0.4) -> su akışını daralt
    model.reactions.EX_h2o_e.bounds = (-1.0, 1.0)

    # Bakım enerjisi: ATPM sabit bir değer (alt=üst) olduğu için iki sınırı BİRLİKTE değiştiriyoruz
    # radyasyon hasarını onarmak ek ATP gerektirir -> varsayılan olarak "bakim_carpani" kat artır
    yeni_bakim = atpm.lower_bound * bakim_carpani
    atpm.bounds = (yeni_bakim, yeni_bakim)
    print(f"Yeni bakım gereksinimi: {atpm.bounds}")

    return model


def mars_buyume(model):
    mars_solution = model.optimize()
    print(f"Mars koşulunda büyüme oranı (1/saat): {mars_solution.objective_value}")
    return mars_solution


def main():
    model = modeli_yukle()
    baseline = referans_buyume(model)
    atpm = bakim_reaksiyonunu_bul(model)
    model = mars_kisitlarini_uygula(model, atpm)
    mars_solution = mars_buyume(model)

    print()
    print("--- Özet ---")
    print(f"Dünya benzeri: {baseline.objective_value}")
    print(f"Mars koşulu:   {mars_solution.objective_value}")


if __name__ == "__main__":
    main()
