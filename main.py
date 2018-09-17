vitamineb12 =["rs2336573","rs1131603","rs3742801","rs2270655","rs12272669","rs34324219","rs7788053","rs602662","rs1801222","rs41281112", "rs1141321","rs652197","rs1801133"]

test =["rs9803031","rs10267"]

#print(variant_genetique_liste(test))

#print(variant_genetique_liste(vitamineb12))

class Recherche_polymorph():

  def __init__(self, poly):
    self.poly= poly
    import os
    import glob
    chemin = os.getcwd()
    self.liste_fichier = glob.glob(chemin+"/*.vcf")

  def variant_genetique_liste2(self):
    
    import re
    patient_gene = []
    
    for f in self.liste_fichier:
        fichier = open(f, "r")
        patient = fichier.read()
        data = patient.split("\n")
        liste = []
        geneti = []
        genetique =[]
        
        liste = [d.split("\t")for d in data]
      
        regex = "##"
        for l in liste :
          if (re.search(regex,l[0])is not None):
            l.clear()

        geneti = [le for le in liste if le != []]
        geneti = geneti[0:len(geneti)-1]  
        regex1 = "\;"
        
        for ge in geneti:
          if (re.search(regex1, ge[2])is not None) :
            ge[2]=ge[2].split(";")
        
        for ge in geneti:
          genetique.append(ge)
          if ge[2]==list(ge[2]):
            i = 1
            lns =[]
            
            while(i<len(ge[2])):
              ni = list(ge)
              lns.append(ni)
              i+=1
            
            for ni in lns:
              ni[2] = ni[2].pop()
              ni[2]= "".join(ni[2])
              genetique.append(ni)
            
            ge[2]= ge[2].pop()
            ge[2]="".join(ge[2])

        
        id_patient = genetique[0][9]
        
        for g in genetique:
          for ln in self.poly:
            if(g[2]==ln):
              resultat = id_patient+":"+ ln+":"+g[9]
              resultat = resultat.split(":")
              resultat = resultat[0:3]
              patient_gene.append(resultat)
                  

    return(patient_gene) 
  
  def tableau_initial(self):
    import pandas
    tab = self.variant_genetique_liste2()
    p = pandas.DataFrame(tab, columns = ["PATIENTS","ID", "GENOTYPE"])
    vf = p.sort_values(by= ["ID","PATIENTS"], ascending = True)
    return vf


  def new_tableau(self):
    import pandas
    liste3 = []
    tableau = self.variant_genetique_liste2()
    for t in tableau:
      liste3.append(t)
      if t[2] == '1/1':
        liste3.append(t)
    p= pandas.DataFrame(liste3, columns = ["PATIENTS", "ID", "GENOTYPE"])
    yx= p.sort_values(by = ["ID","PATIENTS"], ascending = True)
    return yx


  def nb_polymorphisme (self):
    tableau = self.new_tableau()
    return tableau.groupby("ID").size() #tableau["ID"].value_counts()

  def new_donnee(self):
    tableau = self.new_tableau()
    table = self.nb_polymorphisme()
    table1 = [n for n in table]
    tableau2 =[r for r in tableau['ID']]
    tableau3 = sorted(set(tableau2))
    return table1


  def graphique_ID(self):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    #table = variant_genetique_liste2(poly)
    table = self.nb_polymorphisme ()
  
    plt.hist(table, histtype='bar', bins = 'auto', orientation='vertical')
  #return numpy.histogram(table,bins='auto')
    return plt.savefig('plot1.png')  #table["ID"].hist()


  def camembert(self):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    tableau = self.new_tableau()
    table = self.nb_polymorphisme()
    table1 = [n for n in table]
    tableau2 =[r for r in tableau['ID']]
    tableau3 = sorted(set(tableau2))
    plt.figure(figsize = (8, 8))
    plt.pie(table1, labels = tableau3,autopct="%1.2f pourcents")
    plt.legend()
    plt.title("Fréquences allèliques des polymorphismes")
    return plt.savefig('plot3.png')



#sorted(permet de classer les éléments d'une liste dans un ordre croissant en général ou decroissant).
vitamineb12 = Recherche_polymorph(vitamineb12)

tableauinitial = vitamineb12.tableau_initial()

tableauinfo= vitamineb12.new_tableau() 

frequence = vitamineb12.nb_polymorphisme()

graph = vitamineb12.camembert()

print(tableauinitial)

print(tableauinfo)

print(frequence)

print(graph)