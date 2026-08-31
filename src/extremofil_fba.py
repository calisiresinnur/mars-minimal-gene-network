"""
Ekstremofil karşılaştırması: Salinibacter ruber (iMB631) için Mars FBA analizi.

NEDEN BU MODEL: mars_fba.py'nin B. subtilis (iYO844) analizi bir "sıradan"
(mezofilik, düşük tuza adapte) organizma kullanıyor. Kullanıcının sorusu:
gerçekten strese ÖNCEDEN adapte olmuş bir ekstremofil kullansak sonuç
değişir mi? Salinibacter ruber, doygun tuzlu suda (%20-30 tuz) yaşayan,
gerçek bir aşırı halofilik BAKTERİ (arke değil -- hücre biyolojisi
B. subtilis'le daha karşılaştırılabilir) -- yani Mars'ın "düşük su
aktivitesi" eksenine doğal olarak adapte bir organizma.

Model kaynağı: iMB631 (631 gen, 1459 reaksiyon), Ghosh & Mohapatra 2019,
"A genome-scale metabolic network reconstruction of extremely halophilic
bacterium Salinibacter ruber", PLOS ONE, PMC6508672.
SBML: https://doi.org/10.1371/journal.pone.0216336.s002 (S2 File).

ÖNEMLİ FARK -- iYO844 ile doğrudan karşılaştırılamayan noktalar:
1. İsimlendirme: iMB631 ModelSEED tarzı (M_ex00027 gibi cpd ID'leri),
   BiGG tarzı değil (EX_glc__D_e gibi). Exchange yönü de TERS: iYO844'te
   "metabolit_e --> " (negatif akı = alım), iMB631'de " --> metabolit_e"
   (pozitif akı = alım) -- model.medium sözlüğü bu farkı soyutluyor, kodun
   geri kalanında bunu unutmayın.
2. Referans ortam: iYO844 gerçek bir MİNİMAL ortamda (glikoz + inorganik
   tuzlar) büyüyebiliyor (prototrof). iMB631 modelde birden fazla amino
   asit/vitamin dışarıdan verilmezse büyüme SIFIR çıkıyor (tek bir eksik
   besinle açıklanamıyor -- muhtemelen gerçek çoklu oksotrofi ya da
   otomatik rekonstrüksiyonun (ModelSEED tarzı) tamamlanmamış boşlukları).
   Orijinal makale de saf bir minimal ortam değil, pepton+maya özütü
   içeren zengin/tanımsız bir ortam (MGM) kullanmış ve FBA için tam bir
   akı sınırı tablosu yayınlamamış. Bu yüzden REFERANS_ORTAM burada
   bizim seçtiğimiz, açıkça belgelenmiş bir kalibrasyon: glikoz=1.0,
   her bir amino asit/vitamin=0.1 mmol/gDW/h. Bu, makalenin bildirdiği
   0.297/saat referans büyümeye çok yakın bir sonuç veriyor (0.268/saat)
   -- rastgele seçilmedi, bu değere kalibre edildi.

Bu iki fark nedeniyle mutlak büyüme oranlarını veya eşik değerlerini
iYO844 ile BİREBİR karşılaştırmak yerine, her iki organizmanın KENDİ
referansına göre NİSPİ (relative) davranışını (büyüme oranındaki düşüş
yüzdesi, feasibility eşiğinin şiddet ekseni üzerindeki konumu) karşılaştırmak
daha doğru. Bkz. README > Ekstremofil karşılaştırması.
"""

import gzip
import os

import cobra

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_ONBELLEK = os.path.join(PROJE_KOKU, "data", "models", "iMB631.xml.gz")

# Referans (Salinibacter için "Dünya benzeri", yani MODELİN KENDİ doğal/rahat
# koşulu) ortamı -- bkz. yukarıdaki docstring notu.
GLIKOZ_ID = "M_ex00027"
O2_ID = "M_ex00007"
H2O_ID = "M_ex00001"
ATPM_ID = "M_rxn00062"  # "ATP maintenance requirement", iYO844'teki ATPM'nin analogu
BIOMASS_SAHTE_EX_ID = "M_ex11416"  # ModelSEED artefaktı, gerçek değil -- her zaman kapalı

ORGANIK_ID = [
    "M_ex00039",  # L-Lysine
    "M_ex00051",  # L-Arginine
    "M_ex00066",  # L-Phenylalanine
    "M_ex00069",  # L-Tyrosine
    "M_ex00129",  # L-Proline
    "M_ex00322",  # L-Isoleucine
    "M_ex00104",  # Biotin
    "M_ex00220",  # Riboflavin
    "M_ex00393",  # Folate
    "M_ex00423",  # Vitamin B12r
    "M_ex00635",  # Cbl (B12 türevi)
    "M_ex00793",  # Thiamine phosphate
    "M_ex00794",  # Trehaloz (TRHL) -- halofillerde uyumlu çözünen madde
]

REFERANS_GLIKOZ = 1.0
REFERANS_ORGANIK = 0.1
REFERANS_O2 = 3.15  # modelin kendi varsayılan/gap-filling değeri


def modeli_yukle():
    """iMB631 modelini yerel önbellekten yükler (bkz. mars_fba.py'deki aynı
    Windows-Unicode-yol sorununun çözümü -- burada da geçerli)."""
    with gzip.open(MODEL_ONBELLEK, "rt", encoding="utf-8") as f:
        sbml_metni = f.read()
    model = cobra.io.read_sbml_model(sbml_metni)
    print(f"Model yüklendi (iMB631): {len(model.reactions)} reaksiyon, {len(model.genes)} gen")
    return model


def referans_ortami_uygula(model):
    """Modelin kendi doğal/kalibre edilmiş referans ortamını kurar (bkz.
    docstring). Bu, iYO844'ün varsayılan minimal ortamının FONKSİYONEL
    analogu ama SAYISAL olarak birebir karşılaştırılamaz (bkz. yukarı)."""
    med = model.medium
    med[BIOMASS_SAHTE_EX_ID] = 0.0
    med[GLIKOZ_ID] = REFERANS_GLIKOZ
    for rid in ORGANIK_ID:
        if rid in med:
            med[rid] = REFERANS_ORGANIK
    model.medium = med
    return model


def referans_buyume(model):
    referans_ortami_uygula(model)
    baseline = model.optimize()
    print(f"Salinibacter ruber referans büyüme oranı (1/saat): {baseline.objective_value}")
    return baseline


def bakim_reaksiyonunu_bul(model):
    atpm = model.reactions.get_by_id(ATPM_ID)
    print(f"Bakım reaksiyonu: {atpm.id} ({atpm.name}) | mevcut sınırlar: {atpm.bounds}")
    return atpm


def mars_kisitlarini_uygula(model, atpm, o2_ub=REFERANS_O2, glc_ub=REFERANS_GLIKOZ,
                             organik_ub=REFERANS_ORGANIK, h2o_ub=1000.0,
                             bakim_carpani=1, bakim_taban=None, sessiz=False):
    """iYO844'teki mars_kisitlarini_uygula'nın Salinibacter analoğu.

    DİKKAT: burada "lower_bound" değil "upper_bound" kısıtlanıyor çünkü
    iMB631'de exchange yönü ters (pozitif akı = alım) -- bkz. docstring.
    """
    med = model.medium
    med[BIOMASS_SAHTE_EX_ID] = 0.0
    med[O2_ID] = o2_ub
    med[GLIKOZ_ID] = glc_ub
    for rid in ORGANIK_ID:
        if rid in med:
            med[rid] = organik_ub
    med[H2O_ID] = h2o_ub
    model.medium = med

    taban = bakim_taban if bakim_taban is not None else atpm.upper_bound
    yeni_bakim = taban * bakim_carpani
    atpm.bounds = (yeni_bakim, yeni_bakim)
    if not sessiz:
        print(f"Yeni bakım gereksinimi: {atpm.bounds}")

    return model


def mars_buyume(model):
    sol = model.optimize()
    print(f"Mars koşulunda büyüme oranı (1/saat): {sol.objective_value} | durum: {model.solver.status}")
    return sol


def main():
    model = modeli_yukle()
    baseline = referans_buyume(model)
    atpm = bakim_reaksiyonunu_bul(model)
    mars_kisitlarini_uygula(model, atpm, o2_ub=0.3, glc_ub=0.05, organik_ub=0.01, h2o_ub=1.0, bakim_carpani=3)
    mars_solution = mars_buyume(model)

    print()
    print("--- Özet ---")
    print(f"Referans (Salinibacter'in kendi rahat koşulu): {baseline.objective_value}")
    print(f"Mars koşulu:                                   {mars_solution.objective_value}")


if __name__ == "__main__":
    main()
