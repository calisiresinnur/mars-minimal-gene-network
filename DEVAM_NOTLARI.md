# Proje Devam Notları — Başka Bir Sohbetten Devam Etmek İçin

Son güncelleme: 2026-08-31

Bu dosya, bu projede şu ana kadar yapılan her şeyi, alınan kararları ve
sıradaki adımı özetliyor. Yeni bir Claude Code sohbetine bu dosyayı
gösterip "buradan devam et" demen yeterli.

---

## 1) Proje ne, nerede

- **Konu**: Mars yüzey koşullarında bir mikroorganizmanın (genom-ölçekli
  metabolik model + FBA ile) metabolik olarak canlı kalıp kalamayacağını,
  hangi genlerin bu koşullarda esansiyel olduğunu hesaplamalı incelemek.
- **Hedef**: IAC 2026 · IAF/IAA Space Life Sciences Symposium (A1), Paper ID
  114761. Yazar: Esinnur Çalışır, İstanbul Üniversitesi.
- **Ana repo**: https://github.com/calisiresinnur/mars-minimal-gene-network
- **Yerel dizin**: `C:\Users\Ergün\Belgeler\mars-minimal-gene-network`
- Repo `README.md`'si ÇOK DETAYLI — tüm bulgular, kaynaklar, düzeltmeler
  orada da var. Bu dosya sadece hızlı bir özet/harita.

## 2) Şu ana kadar yapılanlar (kronolojik özet)

1. **Ortam bug'ları düzeltildi** (`src/mars_fba.py`): (a) BiGG'in
   http→https yönlendirmesini cobra takip edemiyordu → model artık
   `data/models/iYO844.xml.gz` altında yerelde önbelleğe alınıyor; (b)
   Windows'ta kullanıcı adındaki "ü" karakteri yüzünden libSBML dosya
   yolunu açamıyordu → `.xml.gz` Python'ın gzip modülüyle açılıp SBML
   içeriği libSBML'e **string** olarak veriliyor (dosya yolu hiç
   geçmiyor). **Bu iki sorun her yeni model eklerken tekrar çıkabilir,
   çözüm aynı kalıyor.**
2. **`mars_kalibrasyon.py`** düzeltildi — su kısıtı (`h2o_cap`) daha önce
   sabit kodlanmıştı, kalibrasyon taramasına hiç dahil değildi; artık
   parametrize.
3. **Literatür taraması** yapıldı (Mars atmosfer verileri: %95.54 CO₂,
   %0.13 O₂, %0.03 H₂O, 0.69 kPa — Cortesão ve ark. 2019; radyasyon dozu:
   0.64 mSv/gün — Hassler ve ark. 2014; NGAM referans değerleri; vb.) —
   tüm kaynaklar `README.md` > "Kaynaklar" ve "Karşılaştırmalı genomik"
   bölümlerinde, PMID/DOI ile.
4. **Duyarlılık analizi** (`src/mars_duyarlilik.py`, 357 senaryo) — B.
   subtilis'te (iYO844) her bakım çarpanı için keskin bir feasibility
   "uçurumu" var: eşiğin altında tamamen infeasible, üstünde büyüme
   şiddetle ~doğrusal artıyor.
5. **Tekli gen silme analizi** (`src/mars_gen_silme.py`) —
   ⚠️ **KRİTİK DERS**: İlk çalıştırmada "3 gen (pabB, menC, menD) Mars'ta
   esansiyellikten çıkıyor" diye YANLIŞ bir bulgu raporlandı. Kök neden:
   solver'ın varsayılan feasibility tolerance'ı (1e-7); çok düşük büyüme
   oranlarında (~0.0027/saat) biyokütle denklemindeki kofaktör akışları
   toleransa çok yakın kalıyor ve sahte-feasible sonuç üretiyor. Tolerance
   1e-9'a çekilince (`SOLVER_TOLERANCE` sabiti eklendi) DÜZELTİLMİŞ bulgu:
   **hiçbir gen esansiyellik durumunu değiştirmiyor** (171 esansiyel gen,
   Dünya ve Mars'ın 3 senaryosunda da aynı). **Bundan sonraki her gen
   silme/essentiality analizinde tolerance MUTLAKA 1e-9'a çekilmeli.**
6. **Karşılaştırmalı genomik** (`src/mars_deg_karsilastirma.py`,
   `data/stres_regulon_genleri.csv`) — B. subtilis'in gerçek stres
   regulonlarından (PerR/oksidatif, SigV/hücre zarfı, ResDE/O₂
   kısıtlaması, ISS uçuş DEG'leri) 29 gen modelde bulunup FBA
   sonuçlarıyla çapraz kontrol edildi. Ayrı bir bulgu: `dltABCD` genleri
   modelde HER koşulda esansiyel çıkıyor ama gerçek literatür bu genlerin
   silinmesinin ölümcül OLMADIĞINI gösteriyor (sadece lizozime hassasiyet
   artıyor) — iYO844'ün bilinen bir eksikliği (düz/sübstitüe-edilmemiş
   teikoik asit üreten alternatif reaksiyon yok).
7. **Bilimsel doğruluk denetimi** (kullanıcı isteğiyle) — tüm proje
   literatürle çapraz kontrol edildi. 4 yanlış yazar-adı atfı düzeltildi
   (Moeller→Nicholson, Nickerson→Morrison/Fajardo-Cavazos/Nicholson,
   Newcombe→Cortesão, genel "Nakano" atfı Puri-Taneja 2007 ile
   güçlendirildi). Bu denetim sırasında madde 5'teki kritik hata
   yakalandı ve düzeltildi.
8. **Ekstremofil karşılaştırması** (`src/extremofil_fba.py`,
   `extremofil_duyarlilik.py`, `extremofil_gen_silme.py`) — kullanıcının
   sorusu üzerine: "Dünya'da yaşayan minimal gen seti Mars'ta da işe
   yarıyor" bulgusu organizma seçiminden bağımsız mı? **Salinibacter
   ruber** (iMB631, gerçek aşırı halofilik bakteri, Ghosh & Mohapatra
   2019 PLOS ONE, PMC6508672) eklendi, aynı pipeline paralel çalıştırıldı.
   - Model kalibrasyonu: iMB631'in ham ortamı tüm amino asit/vitaminleri
     sınırsız veriyordu (anlamsız büyüme, 13.88/saat) → elle kalibre
     edildi (glikoz=1.0, her organik kaynak=0.1 mmol/gDW/h) → 0.2676/saat
     (makalenin bildirdiği 0.297/saat'e yakın).
   - **Bulgu**: B. subtilis'in aksine Salinibacter'de test edilen HİÇBİR
     noktada infeasible olmuyor — büyüme tamamen doğrusal, bakım
     çarpanından neredeyse bağımsız (B. subtilis'teki keskin uçurum yok).
   - **Ortak nokta**: yine de gen-esansiyellik seti değişmiyor (148=148,
     her koşulda aynı) — B. subtilis'teki düzeltilmiş sonuçla aynı
     kalitatif bulgu.
   - **Dürüstlük notu**: bu dayanıklılık farkının ne kadarı gerçek
     halofil biyolojisi, ne kadarı iMB631'in daha az olgun (ModelSEED
     tarzı otomatik) rekonstrüksiyon kalitesinden kaynaklanan yapay
     esneklik olduğu bu analizle ayırt edilemiyor — makalede belirtilmeli.

## 3) Önemli metodolojik dersler (tekrar karşılaşabilirsin)

- **Windows Unicode kullanıcı adı** → libSBML dosya yolu okuyamıyor →
  her zaman gzip ile açıp SBML'i **string** olarak `cobra.io.read_sbml_model`'e ver.
- **Bash `/tmp` yolu** Windows Python ile uyuşmuyor (Git Bash path
  translation sorunu) → geçici dosyalar için proje içi göreli yol kullan
  (scratchpad dizini de bir seçenek ama proje-dizini-göreli daha güvenli
  oldu bu oturumda).
- **FBA'da çok düşük büyüme oranlarında solver tolerance kritik** —
  varsayılan (1e-7) sahte-feasible sonuç verebilir. Herhangi bir gen
  silme/essentiality analizinde `model.solver.configuration.tolerances.feasibility = 1e-9`
  satırını EN BAŞTA ekle, özellikle WT büyüme ≲0.01 gibi küçükse.
  Şüpheli durum işareti: KO büyümesi WT büyümesiyle "tesadüfen" neredeyse
  birebir aynıysa (12+ basamak) bu bir tolerance artefaktı olabilir.
- **Yeni bir GEM eklerken** (farklı organizma): (1) exchange reaksiyon
  isimlendirmesi ve YÖNÜ (BiGG: negatif=alım; ModelSEED tarzı: pozitif=
  alım olabilir) kontrol edilmeli — `model.medium` getter/setter'ı bu
  farkı soyutlar, ondan faydalan; (2) modelin "varsayılan" ortamı
  gerçekçi olmayabilir (tümü sınırsız açık) — mutlaka kontrol edip
  gerekirse elle kalibre et; (3) literatürdeki "referans büyüme" değeriyle
  karşılaştırarak kalibrasyonu doğrula.
- **Model seçimi (organizma + kürasyon kalitesi) sonucu kökten
  değiştiriyor** — bu artık ampirik olarak gösterildi (B. subtilis vs
  Salinibacter karşılaştırması).
- Git commit/push için kullanıcı baştan blanket onay verdi ("çalıştır,
  githuba ve bilgisayarıma kaydet") — her adımda tekrar onay istemeye
  gerek yok, ama önemli/riskli bir değişiklik varsa yine de belirt.

## 4) Repo yapısı (mars-minimal-gene-network)

```
.
├── README.md                          # HER ŞEYİN detaylı belgesi — önce buraya bak
├── DEVAM_NOTLARI.md                    # bu dosya
├── requirements.txt                    # cobra, pandas, matplotlib
├── data/
│   ├── models/
│   │   ├── iYO844.xml.gz               # B. subtilis modeli (BiGG)
│   │   └── iMB631.xml.gz               # Salinibacter ruber modeli (PLOS ONE S2 dosyası)
│   └── stres_regulon_genleri.csv       # Literatür stres-regulonu/DEG referans listesi
├── results/                            # tüm CSV/PNG çıktılar
└── src/
    ├── mars_fba.py                     # B. subtilis: model yükleme + Mars kısıtları + FBA
    ├── mars_kalibrasyon.py             # B. subtilis: hızlı senaryo kontrolü
    ├── mars_duyarlilik.py              # B. subtilis: duyarlılık analizi + grafik
    ├── mars_gen_silme.py               # B. subtilis: tekli gen silme (SOLVER_TOLERANCE=1e-9)
    ├── mars_deg_karsilastirma.py       # B. subtilis: FBA x literatür DEG karşılaştırma
    ├── extremofil_fba.py               # Salinibacter: model yükleme + Mars kısıtları + FBA
    ├── extremofil_duyarlilik.py        # Salinibacter: duyarlılık analizi + grafik
    └── extremofil_gen_silme.py         # Salinibacter: tekli gen silme (SOLVER_TOLERANCE=1e-9)
```

## 5) Şu anki durum / kaldığımız yer

Kullanıcı şu soruyu sordu: "Dünya'da yaşayan minimal gen seti Mars'ta da
gerekli" bulgusu B. subtilis'e mi özgü? → Ekstremofil karşılaştırmasıyla
kısmen cevaplandı (madde 2.8). Sonraki adım olarak kullanıcı şunu istedi:

> **Minimal sentetik hücre yaklaşımı** (JCVI-syn3.0/syn3A gibi, ~473-493
> genlik) — "B. subtilis'in hangi genleri gerekli" yerine "Mars'ta
> hayatta kalmak için teorik olarak en az kaç/hangi gen yeterli" sorusunu
> sormak. Kullanıcı bunun kendi asıl ilgilendiği yöne daha yakın
> olduğunu belirtti.

**Kararlar (kullanıcıdan alındı)**:
- Bu **TAMAMEN YENİ, AYRI bir repo/proje** olacak — mevcut
  mars-minimal-gene-network reposuna dahil edilmeyecek. **Henüz
  oluşturulmadı** (ne yerel klasör ne GitHub reposu var).
- Model taramasını ben (Claude) yapacaktım.

**Araştırma bulguları (JCVI-syn3A GEM taraması, tamamlanmadı)**:
- Gerçek, yayınlanmış bir GEM VAR: **Breuer ve ark. 2019, eLife,
  "Essential metabolism for a minimal cell"**
  (https://elifesciences.org/articles/36842, PMC6609329). JCVI-syn3A
  için: **155 gen, 338 reaksiyon, 304 metabolit**. Deneysel transpozon
  mutajenez verisiyle doğrulanmış (in vivo esansiyellik %92, in silico
  %79, Matthews correlation 0.59). Çift. süresi ~2 saat (belirsizlik
  payıyla). **Zengin/tanımsız ortam** kullanıyor -- "tüm biyokütle
  öncüllerini sağlayan zengin bir in silico ortam, glikoz tek enerji
  kaynağı" varsayılmış; gerçek deneysel ortam (SP4 medium) tanımlı/minimal
  değil, yazarlar da "normal büyümeyi destekleyen tanımlı bir ortam henüz
  elde edilmedi" diyor.
  ⚠️ **SBML/model dosyasının tam indirme linki HENÜZ BULUNAMADI** —
  makalenin ek dosyalarına (Supplementary files) veya
  https://github.com/Luthey-Schulten-Lab/Minimal_Cell reposuna
  bakılmalı (bu repo whole-cell kinetic model için, metabolik ağı da
  içeriyor olabilir).
- İkinci bir aday: **iJL208** modeli — "Genome-scale metabolic modeling
  reveals key features of a minimal gene set" (Rees-Garbutt ve ark.(?)
  2021, Molecular Systems Biology, PMC8290834) — JCVI-syn3.0'ı ebeveyni
  JCVI-syn1.0 ile karşılaştıran bir GEM. Bu da incelenmeli, SBML
  kaynağı/indirme linki bulunmalı.
- **Henüz yapılmadı**: hangi modelin (Breuer 2019 metabolik ağı mı,
  iJL208 mi) daha uygun/erişilebilir olduğuna karar verilmedi; SBML
  dosyası indirilip cobra ile test edilmedi.

## 6) Sıradaki somut adımlar — GÜNCELLENDİ 2026-08-31

**Madde 1-4 TAMAMLANDI, ayrı bir projede devam ediyor:**
`C:\Users\Ergün\Belgeler\mars-minimal-cell-network` (henüz GitHub'a
push edilmedi — kullanıcıdan izin bekleniyor/isteniyor). Detaylar için o
projenin kendi `DEVAM_NOTLARI.md`'sine bak. Özet: model = Breuer ve ark.
2019 (eLife 36842) JCVI-syn3A metabolik ağı (iMMSYN, 155 gen/338 rxn/304
metabolit), indirildi ve doğrulandı (eLife'ın "Figures and data"
sayfasından Supplementary file 9); NGAM yapısı (ATPase/Protein_degrad/
RNA_degrad sabit alt sınırları: 0.575/0.00035/0.0077) makalenin tam
metniyle (JATS XML) birebir doğrulandı — ad hoc varsayım gerekmedi;
`src/mars_fba.py` yazıldı ve bir solver-warm-start artefaktı canlı
yakalanıp düzeltildi (infeasible durumda objective_value asla
raporlanmıyor); ilk bulgu: B. subtilis'ten farklı olarak burada asıl
kısıtlayıcı su değil **glikoz** (keskin uçurum glc_lb ≈ -0.8/-0.75
mmol/gDW/h arasında; su tek başına kısıtlayıcı değil, feasible).

**Kalan adımlar (yeni projede, kendi DEVAM_NOTLARI.md'sinin madde 5-6'sı):**

5. Tam duyarlılık analizi + gen esansiyellik/silme analizi
   (SOLVER_TOLERANCE=1e-9'dan başlayarak) — henüz yapılmadı.
6. Üç modelin (B. subtilis/Salinibacter/JCVI-syn3A) karşılaştırmalı
   "kısıtlayıcı darboğaz" bulgusu bir araya getirilmeli.
7. GitHub'a push için kullanıcı onayı (yeni repo oluşturma, mevcut repoya
   commit/push'tan farklı bir eylem olarak ayrıca teyit edilecek).

## 7) Genel hatırlatmalar

- Kullanıcı Türkçe konuşuyor, teknik yazım tarzı: dürüst, aşırı iddialı
  olmayan, kaynaklı, "bulunamadı"yı da rapor eden bir üslup benimsendi —
  bunu koru.
- Her önemli kod değişikliğinden sonra README güncellenip commit/push
  yapılıyor (kullanıcı blanket onay verdi).
- Kullanıcı **bilimsel doğruluğa çok önem veriyor** — literatür
  atıflarını ASLA hafızadan yazma, her zaman WebSearch/WebFetch ile
  çapraz doğrula (bu oturumda 4 kez yanlış yazar adı hatırlanıp
  düzeltildi).
