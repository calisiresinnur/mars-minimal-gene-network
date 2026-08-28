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
└── src/
    └── mars_fba.py      # Model yükleme + Mars kısıtları + FBA (ana script)
```

## Kurulum ve çalıştırma

Bu adımları VS Code'un içindeki terminalde (Ctrl+`) sırayla çalıştır:

```
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
python src/mars_fba.py
```

PowerShell "çalıştırma politikası" hatası verirse (venv aktive olmuyorsa), önce şunu
çalıştır ve tekrar dene:

```
Set-ExecutionPolicy -Scope CurrentUser RemoteSigned
```

## Durum

- [x] Referans (Dünya benzeri) büyüme doğrulandı — 0.118 /saat
- [x] Mars kısıtları tanımlandı (O₂, CO₂, organik karbon, su, bakım enerjisi)
- [ ] Tekli gen silme (single gene deletion) analizi
- [ ] Karşılaştırmalı genomik (DEG, stres-toleransı gen listeleri) entegrasyonu
- [ ] Sonuç şekilleri ve tablolar
- [ ] Tam metin yazımı

Ayrıntılı 14 günlük yol haritası ve kaynak linkleri için: proje sohbetindeki
Mars Gen Ağı Yol Haritası dokümanına bakabilirsin.
