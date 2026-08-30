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

Model kaynağı: iYO844, BiGG Models (http://bigg.ucsd.edu/models/iYO844).
cobra.io.load_model("iYO844") BiGG'nin canlı web API'sine bağımlıdır ve BiGG
http->https yönlendirmesi kurulu cobra/httpx sürümüyle başarısız olabiliyor
(bkz. modeli_yukle). Bu yüzden model dosyasını data/models/ altında yerelde
önbelleğe alıyoruz; script hem çevrimdışı hem de tekrarlanabilir çalışır.

Ayrıca: Windows'ta kullanıcı adında ASCII-dışı karakter varsa (ör. "Ergün")
libSBML'in dosya yolunu C seviyesinde açması başarısız oluyor ("No SBML model
detected in file" hatası) — path'in kendisi geçerli olsa da. Bunu aşmak için
.xml.gz dosyasını libSBML'e hiç yol vermeden, Python'ın kendi gzip modülüyle
açıp ham SBML metnini string olarak cobra'ya veriyoruz.
"""

import gzip
import os
import urllib.request

import cobra

PROJE_KOKU = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_ONBELLEK = os.path.join(PROJE_KOKU, "data", "models", "iYO844.xml.gz")
MODEL_URL = "https://bigg.ucsd.edu/static/models/iYO844.xml.gz"


def modeli_yukle():
    """iYO844 modelini yükler.

    Öncelik sırası:
      1) Yerel önbellek (data/models/iYO844.xml.gz) varsa doğrudan oradan oku.
      2) Yoksa BiGG'den https üzerinden indirip önbelleğe kaydet, sonra oku.
      3) O da başarısız olursa cobra'nın kendi çevrimiçi yükleyicisini dene
         (BiGG'in http->https yönlendirmesini takip edemeyen eski sürümlerde
         bu adım başarısız olabilir; hata mesajı bunu açıkça belirtir).
    """
    if not os.path.exists(MODEL_ONBELLEK):
        os.makedirs(os.path.dirname(MODEL_ONBELLEK), exist_ok=True)
        try:
            print(f"Model önbellekte yok, indiriliyor: {MODEL_URL}")
            urllib.request.urlretrieve(MODEL_URL, MODEL_ONBELLEK)
        except Exception as indirme_hatasi:
            print(f"Doğrudan indirme başarısız ({indirme_hatasi}), cobra.io.load_model deneniyor...")
            try:
                model = cobra.io.load_model("iYO844")
                print(f"Model yüklendi: {len(model.reactions)} reaksiyon, {len(model.genes)} gen")
                return model
            except Exception as cobra_hatasi:
                raise RuntimeError(
                    "iYO844 modeli ne yerel önbellekten ne de BiGG'den yüklenebildi. "
                    f"İndirme hatası: {indirme_hatasi} | cobra.io.load_model hatası: {cobra_hatasi}"
                ) from cobra_hatasi

    # Not: dosya yolunu doğrudan cobra.io.read_sbml_model'e vermiyoruz -- libSBML,
    # yolda ASCII-dışı karakter (ör. "Ergün") olduğunda dosyayı açamıyor. Bunun
    # yerine gzip'i Python'da açıp ham SBML metnini string olarak veriyoruz.
    with gzip.open(MODEL_ONBELLEK, "rt", encoding="utf-8") as f:
        sbml_metni = f.read()
    model = cobra.io.read_sbml_model(sbml_metni)
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


def mars_kisitlarini_uygula(
    model, atpm, o2_lb=-0.5, glc_lb=-0.05, h2o_cap=1.0, bakim_carpani=3, bakim_taban=None, sessiz=False
):
    # O2: Mars atmosferinin sadece ~%0.17'si O2 + toplam basınç Dünya'nın ~%0.6'sı -> ciddi kısıtla
    model.reactions.EX_o2_e.lower_bound = o2_lb

    # CO2: Mars atmosferinin ~%95'i CO2 -> alıma açıkça izin ver, kısıtlayıcı olmasın
    model.reactions.EX_co2_e.bounds = (-1000, 1000)

    # Organik karbon: Mars yüzeyinde serbest glikoz/organik karbon yok denecek kadar az
    model.reactions.EX_glc__D_e.lower_bound = glc_lb

    # Su: düşük su aktivitesi (aw ~0.4) -> su akışını daralt.
    # DİKKAT: bu sınır metabolik su ÜRETİMİNİ de (solunumun yan ürünü) kapsıyor, sadece
    # dışarıdan su alımını değil. Referans (Dünya benzeri) optimumda hücre büyüme
    # OLMADAN, sadece sabit bakım gereksinimini (ATPM) karşılamak için bile +7.87
    # birim su atımı yapıyor -- yani h2o_cap=1.0 tek başına, O2/glikoz ne olursa olsun
    # modeli infeasible yapan bir kısıt (bkz. kalibrasyon script'indeki bulgular).
    model.reactions.EX_h2o_e.bounds = (-h2o_cap, h2o_cap)

    # Bakım enerjisi: ATPM sabit bir değer (alt=üst) olduğu için iki sınırı BİRLİKTE değiştiriyoruz
    # radyasyon hasarını onarmak ek ATP gerektirir -> varsayılan olarak "bakim_carpani" kat artır.
    #
    # DİKKAT (gerçek bir üretim hatasından öğrenildi): çarpanı ATPM'nin O ANKİ
    # lower_bound'una uygularsak, aynı model/atpm nesnesi birden çok kez (ör. bir
    # tarama döngüsünde) bu fonksiyondan geçirildiğinde değer her seferinde
    # katlanarak büyür (9 -> 13.5 -> 20.25 -> ... sonsuza). Bunu önlemek için
    # taban değeri HER ZAMAN açıkça bilinen bir referanstan alıyoruz: çağıran
    # "bakim_taban" vermezse, atpm.lower_bound'u SADECE bu fonksiyon hiç
    # çağrılmamışsa (yani hâlâ modelin orijinal değeriyse) taban olarak kabul
    # ediyoruz -- modeli tekrar kullanan çağıranlar (ör. mars_duyarlilik.py)
    # bakim_taban'ı MUTLAKA açıkça vermeli.
    taban = bakim_taban if bakim_taban is not None else atpm.lower_bound
    yeni_bakim = taban * bakim_carpani
    atpm.bounds = (yeni_bakim, yeni_bakim)
    if not sessiz:
        print(f"Yeni bakım gereksinimi: {atpm.bounds}")

    return model


def mars_buyume(model):
    mars_solution = model.optimize()
    durum = model.solver.status
    print(f"Mars koşulunda büyüme oranı (1/saat): {mars_solution.objective_value} | durum: {durum}")
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
