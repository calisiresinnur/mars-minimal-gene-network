# Mars Yüzey Koşulları İçin Minimal Gen Ağı Modellemesi

**IAC 2026 · IAF/IAA Space Life Sciences Symposium (A1) · Paper ID 114761**
Yazar: Esinnur Çalışır, İstanbul Üniversitesi

## Proje ne yapıyor

Bu proje, kürasyonu yapılmış bir genom-ölçekli metabolik model (GEM) üzerinden, Mars
yüzey koşullarını (düşük atmosfer basıncı, kısıtlı O₂, bol CO₂, düşük su aktivitesi,
yüksek radyasyon) sayısal kısıtlara çevirip Flux Balance Analysis (FBA) ile
"bu mikroorganizma Mars'ta metabolik olarak canlı kalabilir mi, hangi genler bu
koşullarda esansiyel hale geliyor?" sorusunu hesaplamalı olarak inceliyor.

Metabolik iskelet: **iYO844** (*Bacillus subtilis* genom-ölçekli metabolik modeli,
BiGG/BioModels üzerinden hazır ve kürasyonu yapılmış).

## Klasör yapısı

```
.
├── README.md
├── requirements.txt
├── .gitignore
├── data/
│   └── models/
│       └── iYO844.xml.gz          # Yerel model önbelleği (bkz. "Model kaynağı")
├── results/
│   ├── duyarlilik_sonuclari.csv   # mars_duyarlilik.py çıktısı (357 senaryo)
│   └── buyume_vs_siddet.png       # Büyüme oranı / kısıt şiddeti grafiği
└── src/
    ├── mars_fba.py           # Model yükleme + Mars kısıtları + FBA (ana script)
    ├── mars_kalibrasyon.py   # İlk elle-seçilmiş senaryolarla hızlı doğruluk kontrolü
    └── mars_duyarlilik.py    # Sistematik duyarlılık analizi + grafik/CSV çıktısı
```

## Kurulum ve çalıştırma

Bu adımları VS Code'un içindeki terminalde (Ctrl+`) sırayla çalıştır:

```
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
python src/mars_fba.py
python src/mars_kalibrasyon.py
python src/mars_duyarlilik.py
```

PowerShell "çalıştırma politikası" hatası verirse (venv aktive olmuyorsa), önce şunu
çalıştır ve tekrar dene:

```
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

## Model kaynağı ve bilinen ortam sorunları

Metabolik iskelet **iYO844**, BiGG Models'ten (`bigg.ucsd.edu`) alınır ve
`data/models/iYO844.xml.gz` altında yerelde önbelleğe alınır (dosya repoya
commit'lenmiştir, ~290 KB). `modeli_yukle()` önce bu önbelleğe bakar; yoksa
BiGG'den https üzerinden indirip önbelleğe kaydeder.

Bu önbellekleme iki ayrı ortam sorununu aştığı için gerekli:

1. **BiGG http→https yönlendirmesi**: `cobra.io.load_model("iYO844")` BiGG'nin
   canlı API'sini `http://` ile çağırıyor; BiGG artık `https://`'ye 301
   yönlendirme yapıyor ve kurulu cobra/httpx sürümü bu yönlendirmeyi takip
   etmiyor, bağlantı hatasıyla düşüyor.
2. **Windows'ta Unicode kullanıcı adı**: Kullanıcı adında ASCII-dışı karakter
   varsa (ör. `Ergün`), libSBML dosya yolunu C seviyesinde açamıyor ve
   "No SBML model detected in file" hatası veriyor — dosyanın kendisi ve yolu
   tamamen geçerli olsa bile. Çözüm: `.xml.gz` dosyası Python'ın kendi `gzip`
   modülüyle açılıp ham SBML metni libSBML'e **string** olarak veriliyor
   (dosya yolu hiç libSBML'e geçmiyor).

## Kalibrasyon bulguları (ilk çalıştırma, 2026-08-30)

- Referans (Dünya benzeri) büyüme doğrulandı: **0.1180 /saat**.
- İlk Mars varsayımı (O₂ lb=-0.5, glikoz lb=-0.05, su kısıtı ±1, bakım×3)
  **infeasible** (hayat matematiksel olarak imkânsız), sadece "düşük büyüme"
  değil.
- Kök neden: ATPM (bakım enerjisi) reaksiyonu modelde `lower_bound=upper_bound`
  ile **sabit** bir zorunlu akı (varsayılan 9 mmol ATP/gDW/h). Referans
  optimumda hücre, büyüme sıfır olsa bile sırf bu sabit bakımı karşılamak için
  +7.87 birim su atıyor, -5.7 O₂ ve -1.7 glikoz tüketiyor. Bu yüzden O₂, glikoz
  veya su tek başına belli bir eşiğin altına çekildiğinde model tamamen
  infeasible oluyor — kademeli "yumuşak" bir düşüş değil, bir uçurum.
- Tek tek (diğerleri varsayılanken) bulunan eşikler:
  - O₂ alt sınırı: **lb=-1.0 infeasible → lb=-2.0 feasible** (büyüme≈0.025)
  - Glikoz alt sınırı: **lb=-0.2 infeasible → lb=-0.5 feasible** (büyüme≈0.011)
  - Su kısıtı: **cap=±1 infeasible → cap=±2 feasible** (büyüme≈0.016)
  - Bakım çarpanı tek başına: ×3'e kadar sorun çıkarmıyor (büyüme≈0.051)
- Önemli: su kısıtı önceki kod sürümünde kalibrasyon taramasına hiç dahil
  edilmemiş, sabit ±1 olarak kalmıştı — bu da ilk 5 senaryonun tamamının
  infeasible çıkmasının asıl (ve tek taranmamış) nedeniydi.
- Eşikler ayrı ayrı bulunsa da kısıtlar **etkileşimli**: tek başına feasible
  olan eşik değerleri bir arada uygulandığında (ör. O₂=-2, glikoz=-0.5,
  su=±2, bakım×1.5) yine infeasible çıkabiliyor.

## Kaynaklar (literatür taraması, 2026-08-30)

Kısıt değerlerini gerekçelendirmek için yapılan taramanın bulguları:

- **Gerçek Mars atmosfer bileşimi**: %95.54 CO₂, %0.13 O₂, %0.03 H₂O buharı,
  toplam basınç 0.69 kPa — ölçülmüş değerler.
  [Bacillus subtilis Spore Resistance to Simulated Mars Surface Conditions](https://pmc.ncbi.nlm.nih.gov/articles/PMC6399134/)
  (*Frontiers in Microbiology*, PMC6399134). Aynı çalışmadan: UV'den korunan
  sporlarda canlılık ~%73 korunmuş, UV dahil edildiğinde (8 saatte 115 kJ/m²
  UV-C) ~%6.6'ya düşmüş — Mars'ta asıl öldürücü faktörün kozmik radyasyondan
  çok **UV** olduğuna işaret ediyor.
- **Kozmik (iyonlaştırıcı) radyasyon dozu**: NASA Curiosity/MSL-RAD ölçümleri,
  yüzeyde günlük ~0.64-0.67 mSv, yıllık ~150-250 mSv (Dünya arka planının,
  ~2-3 mSv/yıl, üstünde ama kronik bir doz — akut hücre ölümüne yol açacak
  düzeyde değil).
- **NGAM/ATPM referans değerleri**: BioNumbers'ta tür-spesifik "normal" bakım
  enerjisi değerleri var (ör. *Geobacter metallireducens* ~0.81 mmol
  ATP/gDW/h) ama radyasyon altında bu değerin nasıl arttığına dair
  sayısallaştırılmış bir kaynak **yok**.
- **Sonuç**: [Genome-scale metabolic modelling of extremophiles and its applications in astrobiological environments](https://pmc.ncbi.nlm.nih.gov/articles/PMC10866088/)
  (Noirungsee ve ark. 2024, *Environ Microbiol Reports*) gibi güncel bir
  derleme bile "radyasyon etkileri analize dahil edilebilir" diyor ama somut
  bir sayısal dönüştürme yöntemi önermiyor — yani **radyasyon → ATP bakım
  maliyeti çarpanı için literatürde gerekçelendirilebilir tek bir sayı yok.**
  Bu, alanın kendisinin henüz kapatmadığı bir boşluk.

**Metodolojik karar (bu boşluk nedeniyle):** Tek bir "kanonik" Mars senaryosu
iddia etmek yerine kısıt şiddeti bir **duyarlılık analizi (sensitivity
analysis)** olarak ele alınıyor — bkz. `mars_duyarlilik.py` ve aşağıdaki
grafik. O₂/CO₂/su yüzdeleri yukarıdaki ölçümlerle *niteliksel* olarak
gerekçelendirilebilir (Mars'ta bunlar gerçekten çok kısıtlı) ama bu
yüzdeleri doğrudan bir FBA akı sınırına (mmol/gDW/h) çevirecek bir
kinetik/taşınım modeli yok — böyle bir dönüşüm yapmak sahte bir kesinlik
iddiası olurdu. Bakım çarpanı için de aynı durum geçerli. Bu yüzden her ikisi
de makalede **taranan parametreler** olarak sunulmalı, tek bir "doğru" değer
olarak değil.

## Duyarlılık analizi (`mars_duyarlilik.py`)

O₂, glikoz ve su kısıtları ortak bir şiddet ekseninde (`t`: 0 = en sert ilk
varsayım, 1 = ılımlı bir uç) birlikte taranıyor; bu, 7 farklı bakım
çarpanı (×1.0 – ×4.0) için tekrarlanıyor (357 senaryo, `results/duyarlilik_sonuclari.csv`).

![Büyüme oranı / kısıt şiddeti](results/buyume_vs_siddet.png)

Gözlem: her bakım çarpanı için belirli bir şiddet eşiğine kadar model
**tamamen infeasible** (bir uçurum), eşiğin hemen ötesinde büyüme oranı
şiddetle yaklaşık doğrusal artıyor. Bakım çarpanı arttıkça hem eşik sağa
kayıyor (hayatta kalmak için daha fazla kaynak gerekiyor) hem de aynı şiddet
seviyesindeki maksimum büyüme oranı düşüyor — beklenen, tutarlı bir davranış.

## Durum

- [x] Referans (Dünya benzeri) büyüme doğrulandı — 0.118 /saat
- [x] Mars kısıtları tanımlandı (O₂, CO₂, organik karbon, su, bakım enerjisi)
- [x] Model yükleme ortam hataları giderildi (BiGG yönlendirme + Windows Unicode yol sorunu)
- [x] Kısıt şiddeti kalibrasyonu — ilk feasible Mars senaryosu bulundu
- [x] Literatür taraması — kısıt değerlerinin gerekçelendirilebilirlik sınırları belgelendi
- [x] Sistematik duyarlılık analizi + sonuç grafiği (`mars_duyarlilik.py`)
- [ ] Tekli gen silme (single gene deletion) analizi
- [ ] Karşılaştırmalı genomik (DEG, stres-toleransı gen listeleri) entegrasyonu
- [ ] Tam metin yazımı

Ayrıntılı 14 günlük yol haritası ve kaynak linkleri için: proje sohbetindeki
Mars Gen Ağı Yol Haritası dokümanına bakabilirsin.
