"""
Mars kısıt şiddeti kalibrasyonu — kademeli senaryo taraması.

İlk varsayım seti (O2 -0.5, glikoz -0.05, bakım x3) modeli 'infeasible' yaptı —
yani hayat tamamen imkânsız hale geldi. Bu script kısıtları kademeli gevşeterek
"hangi şiddete kadar hayat hâlâ mümkün" sorusunu haritalıyor.

Her senaryo için modeli sıfırdan yüklüyoruz ki önceki denemenin izleri karışmasın.

Model yükleme ve kısıt uygulama mantığı mars_fba.py ile ortak; burada tekrar
tanımlamak yerine oradan import ediyoruz.

TEŞHİS NOTU (ilk çalıştırmadan sonra eklendi):
Önceki sürümde su kısıtı (h2o_cap) sabit ±1.0 olarak koda gömülüydü ve hiçbir
senaryoda değişmiyordu. Tek başına test edildiğinde ±1.0 su kısıtı, O2 ve glikoz
tamamen serbest bırakılsa BİLE modeli infeasible yapıyor -- çünkü referans
(Dünya benzeri) optimumda hücre, büyüme sıfır olsa bile sadece sabit bakım
gereksinimini (ATPM=9) karşılamak için +7.87 birim su atmak zorunda. Yani ilk 5
senaryonun tamamının infeasible çıkmasının asıl nedeni O2/glikoz/bakım değil,
hiç taranmayan bu su kısıtıydı. Bu script artık su kısıtını da 4. parametre
olarak tarıyor ve tek-tek eşik taramasıyla bulunan gerçek sınırları içeriyor:
  - O2:     lb=-1.0 infeasible  | lb=-2.0 optimal (büyüme≈0.025)
  - Glikoz: lb=-0.2 infeasible  | lb=-0.5 optimal (büyüme≈0.011)
  - Su:     cap=±1  infeasible  | cap=±2  optimal (büyüme≈0.016)
  (diğer üç kısıt varsayılanında sabitken, tek tek taranarak bulundu)
"""

from mars_fba import bakim_reaksiyonunu_bul, mars_kisitlarini_uygula, modeli_yukle


def senaryo_calistir(o2_lb, glc_lb, h2o_cap, bakim_carpani, etiket):
    model = modeli_yukle()  # HER SENARYO İÇİN TEMİZ BİR MODEL
    atpm = bakim_reaksiyonunu_bul(model)
    mars_kisitlarini_uygula(
        model, atpm, o2_lb=o2_lb, glc_lb=glc_lb, h2o_cap=h2o_cap, bakim_carpani=bakim_carpani
    )

    sol = model.optimize()
    durum = model.solver.status
    print(f"{etiket:38s} | durum: {durum:12s} | büyüme: {sol.objective_value}")
    return durum, sol.objective_value


def main():
    print("Senaryo taraması — kademeli gevşetme (artık su kısıtı da dahil)\n")
    senaryo_calistir(-0.5, -0.05, 1.0, 3, "1) İlk varsayım (çok sert)")
    senaryo_calistir(-2.0, -0.05, 1.0, 3, "2) O2 gevşetildi")
    senaryo_calistir(-2.0, -0.2, 1.0, 3, "3) + organik karbon gevşetildi")
    senaryo_calistir(-2.0, -0.2, 1.0, 1.5, "4) + bakım enerjisi gevşetildi")
    senaryo_calistir(-5.0, -0.5, 1.0, 1.5, "5) Daha ılımlı Mars senaryosu (eski su kısıtıyla)")
    senaryo_calistir(-5.0, -0.5, 2.0, 1.5, "6) Aynısı, su kısıtı eşiğinde (±2)")
    senaryo_calistir(-2.0, -0.5, 2.0, 1.5, "7) Eşik değerlerinde sert Mars senaryosu")
    senaryo_calistir(-5.0, -1.0, 8.0, 2.0, "8) Ilımlı-orta Mars senaryosu")


if __name__ == "__main__":
    main()
