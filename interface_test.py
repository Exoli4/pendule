import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import os
import shutil
from math import ceil
import tkinter.font as tkFont
import ctypes
ctypes.windll.shcore.SetProcessDpiAwareness(1)

g=10

def classement(liste):
    return [sorted(liste).index(x) + 1 for x in liste]
def Matrice(N,M,données,limites=0):   #N le nombre de rangées, M le nombre de colonnes
    diag, diag1 = données[:,1],données[:,0] #diag:onnées sur la diagonale principale, 
                                            #diag1:données juste en dessous ou juste au dessus de la diagonale supérieure
    k,T1=0,0 #indice pour la valeur de diaga ou diag1

    mat=[]
    for i in range(N):
        rangée_n=[]
        for j in range(M):
            if i==j:
                rangée_n.append((diag[k]))
                k+=1  #alterne si plus d'une valeur
            elif i==j+1 :
                rangée_n.append((diag1[T1]))
                T1+=1 #alterne si plus d'une valeur
            elif i==j-1:
                rangée_n.append((diag1[T1]))               
            else:
                rangée_n.append(0)

            #conditions limites
            if (j==i and j==0):
                rangée_n[j]+=limites[0]
            if  (j==i and j==N-1):
                if N%2==0:
                    rangée_n[j]+=limites[0]
                else:
                    rangée_n[j]+=limites[1]
            if (i==0 and j==N-1) or (j==0 and i==N-1):
                rangée_n[j]+=limites[2]

            if k+1>len(diag):
                k=0    #réinitialise la valeur pour alterner
            if T1+1>len(diag1):
                T1=0   #réinitialise la valeur pour alterner      
        mat.append(rangée_n)
    return mat
def v_propres(mat,données,limites,N):
    m=mat(N,N,données,limites)
    w,amp=(np.linalg.eig(np.array(m)))

    for i in range((w).size):
        if np.abs(w[i])<=1e-9:
            w[i]=0
  

    return w,amp
def Calc_amp_modes_f_ext(Mat_vect_propres,w,f_ext_amp,wd,gamma):
    
    fn=np.linalg.inv(Mat_vect_propres)@f_ext_amp
    A=np.zeros((len(w),len(wd)))
    for i in range(len(A)):
        A[i]=np.abs(fn[i])/np.sqrt((w[i]**2-wd**2)**2+(gamma*wd)**2)

    return A

class Menu(tk.Frame):
    def __init__(self, Main, centre):
        super().__init__(Main)
        self.centre = centre
        self.centre.title("Oscillateurs couplés")
        
        # Créer des boutons et centrer les éléments
        label = tk.Label(self, text="Observations sur les oscillateurs couplés", font=("Helvetica", 20))
        label.grid(row=0, column=0, columnspan=3, pady=14, sticky="nsew")

        bouton_simulation = tk.Button(self, text="Pendules couplés (sans force externe)", 
                                      command=lambda: centre.show_page(Pendule), font=("Helvetica", 14))
        bouton_simulation.grid(row=1, column=1, columnspan=1, pady=10, sticky="nsew")
        bouton_simulation['borderwidth'] = 3
        bouton_simulation['relief'] = 'raised'

        bouton_simulation = tk.Button(self, text="Pendules couplés (force externe)", 
                                      command=lambda: centre.show_page(Pendule_fext), font=("Helvetica", 14))
        bouton_simulation.grid(row=2, column=1, columnspan=1, pady=10, sticky="nsew")
        bouton_simulation['borderwidth'] = 3
        bouton_simulation['relief'] = 'raised'
        '''
        bouton_autre_page = tk.Button(self, text="Autre Page", 
                                      command=lambda: centre.show_page(Comparaison_T2))
        bouton_autre_page.grid(row=2, column=1, columnspan=1, pady=10, sticky="nsew")
        bouton_autre_page['borderwidth'] = 3
        bouton_autre_page['relief'] = 'raised'
        '''
        # Configurer les lignes et colonnes de la grille pour que les widgets soient centrés
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)


class Pendule(tk.Frame):
    def __init__(self, Main, centre):
        super().__init__(Main)
        self.centre = centre

        # Créer le contenu de la page simulation
        label = tk.Label(self, text="Pendules couplés (sans force externe)", font=("Helvetica", 20))
        label.grid(row=0, column=1, pady=14, sticky="nsew")

        bouton_retour = tk.Button(self, text="Retour au Menu", 
                                  command=lambda: centre.show_page(Menu), font=("Helvetica", 14))
        bouton_retour.grid(row=0, column=0, pady=10, sticky="nsew")

        
        bouton_T2 = tk.Button(self, text="Comparaison T2 (pas de force externe)", 
                                      command=lambda: centre.show_page(Comparaison_T2), font=("Helvetica", 14))
        bouton_T2.grid(row=1, column=1, columnspan=1, pady=10, sticky="nsew")
        bouton_T2['borderwidth'] = 3
        bouton_T2['relief'] = 'raised'

        bouton_T2 = tk.Button(self, text="Comparaison T2 sans extrémités (pas de force externe)", 
                                      command=lambda: centre.show_page(Comparaison_T2_cercle), font=("Helvetica", 14))
        bouton_T2.grid(row=2, column=1, columnspan=1, pady=10, sticky="nsew")
        bouton_T2['borderwidth'] = 3
        bouton_T2['relief'] = 'raised'

        bouton_T2 = tk.Button(self, text="Comparaison T2 2D (pas de force externe)", 
                                      command=lambda: centre.show_page(Comparaison_T2_2D), font=("Helvetica", 14))
        bouton_T2.grid(row=3, column=1, columnspan=1, pady=10, sticky="nsew")
        bouton_T2['borderwidth'] = 3
        bouton_T2['relief'] = 'raised'
        

        # Configurer les lignes et colonnes pour centrer le contenu
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=2)
        self.grid_rowconfigure(2, weight=2)
        self.grid_rowconfigure(3, weight=2)

        self.grid_columnconfigure(0, weight=2)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=2)

class Pendule_fext(tk.Frame):
    def __init__(self, Main, centre):
        super().__init__(Main)
        self.centre = centre

        # Créer le contenu de la page simulation
        label = tk.Label(self, text="Pendules couplés (avec force externe)", font=("Helvetica", 20))
        label.grid(row=0, column=1, pady=14, sticky="nsew")

        bouton_retour = tk.Button(self, text="Retour au Menu", 
                                  command=lambda: centre.show_page(Menu), font=("Helvetica", 14))
        bouton_retour.grid(row=0, column=0, pady=10, sticky="nsew")

        
        bouton_T2 = tk.Button(self, text="Comparaison T2 (pas de force externe)", 
                                      command=lambda: centre.show_page(Comparaison_f_ext), font=("Helvetica", 14))
        bouton_T2.grid(row=1, column=1, columnspan=1, pady=10, sticky="nsew")
        bouton_T2['borderwidth'] = 3
        bouton_T2['relief'] = 'raised'

        bouton_T2 = tk.Button(self, text="Comparaison T2 sans extrémités (pas de force externe)", 
                                      command=lambda: centre.show_page(Comparaison_T2_cercle), font=("Helvetica", 14))
        #bouton_T2.grid(row=2, column=1, columnspan=1, pady=10, sticky="nsew")
        bouton_T2['borderwidth'] = 3
        bouton_T2['relief'] = 'raised'

        bouton_T2 = tk.Button(self, text="Comparaison T2 2D (pas de force externe)", 
                                      command=lambda: centre.show_page(Comparaison_T2_2D), font=("Helvetica", 14))
        #bouton_T2.grid(row=3, column=1, columnspan=1, pady=10, sticky="nsew")
        bouton_T2['borderwidth'] = 3
        bouton_T2['relief'] = 'raised'
        

        # Configurer les lignes et colonnes pour centrer le contenu
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=2)
        self.grid_rowconfigure(2, weight=2)
        self.grid_rowconfigure(3, weight=2)

        self.grid_columnconfigure(0, weight=2)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=2)


class Comparaison_T2(tk.Frame):
    def __init__(self, Main, centre):
        super().__init__(Main)
        self.centre = centre

        # Créer le contenu de l'autre page
        label = tk.Label(self, text="Comparaison T2 (pas de force externe)", font=("Helvetica", 20))
        label.grid(row=0, column=1, pady=14, sticky="nsew")

        bouton_retour = tk.Button(self, text="Retour en arrière", 
                                  command=lambda: centre.show_page(Pendule), font=("Helvetica", 14))
        bouton_retour.grid(row=0, column=0, pady=10, sticky="nsew")

        # Fonction pour récupérer les valeurs et fermer la fenêtre
        def valider():
            try:
                self.centre.Valeurs['N'] = int(self.entry_sites.get())
                self.centre.Valeurs['longueur'] = float(self.entry_longueur.get())
                #self.centre.Valeurs['T1'] = float(self.entry_T1.get())
                self.centre.Valeurs['masse'] = float(self.entry_masse.get())
                centre.show_page(Graph_comp_T2)
            
            except ValueError:
                messagebox.showerror("Erreur", "Veuillez entrer des valeurs numériques valides.")

        # Création des widgets
        self.label_sites = tk.Label(self, text="Nombres de sites :")
        self.entry_sites = tk.Entry(self)
        self.label_longueur = tk.Label(self, text="Longueur (m):")
        self.entry_longueur = tk.Entry(self)
        self.entry_longueur.insert(0,'1.0')
        #self.label_T1 = tk.Label(self, text="T1 (N/m):")
        #self.entry_T1 = tk.Entry(self)
        #self.entry_T1.insert(0,'15.0')
        self.label_masse = tk.Label(self, text="Masse (kg):")
        self.entry_masse = tk.Entry(self)
        self.entry_masse.insert(0,'1.0')
        bouton_valider = tk.Button(self, text="Valider", command=valider)

        # Placement des widgets
        self.label_sites.grid(row=1, column=0, padx=10, pady=5)
        self.entry_sites.grid(row=1, column=1, padx=10, pady=5)
        self.label_longueur.grid(row=2, column=0, padx=10, pady=5)
        self.entry_longueur.grid(row=2, column=1, padx=10, pady=5)
        #self.label_T1.grid(row=3, column=0, padx=10, pady=5)
        #self.entry_T1.grid(row=3, column=1, padx=10, pady=5)
        self.label_masse.grid(row=3, column=0, padx=10, pady=5)
        self.entry_masse.grid(row=3, column=1, padx=10, pady=5)
        bouton_valider.grid(row=4,column=1, pady=10)

        # Configurer les lignes et colonnes pour centrer le contenu
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=2)
        self.grid_rowconfigure(2, weight=2)
        self.grid_rowconfigure(3, weight=2)
        self.grid_rowconfigure(4, weight=2)

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)

class Graph_comp_T2(tk.Frame):
    def __init__(self, Main, centre):
        super().__init__(Main)
        self.update_after_id = None
        self.centre = centre

        bouton_retour = tk.Button(self, text="Retour en arrière", 
                                  command=lambda: centre.show_page(Comparaison_T2))
        bouton_retour.grid(row=0, column=0, pady=10, sticky="nsew")

        def plus():
            N = self.centre.Valeurs.get('N', 2)
            N+=1
            self.centre.Valeurs['N']=N
            if N<12:
                self.afficher_subplots()
            self.update_plot(None)
            
        def moins():
            N = self.centre.Valeurs.get('N', 2)
            N-=1
            self.centre.Valeurs['N']=N
            if N<12:
                self.afficher_subplots()
            self.update_plot(None)

        def enregistrer():
            N = self.centre.Valeurs.get('N', 2)
            T1 = self.centre.Valeurs.get('T1', 2)
            T2 = self.centre.Valeurs.get('T2', 2)
            chemin = f'C:/Users/Mathieu/.anaconda/STAGE/modes_normaux_1D/{N}'
            if os.path.exists(chemin)==False:
                os.makedirs(chemin)
            chemin_image = os.path.join(chemin, f"T1={T1},T2={T2}.png")    
            self.fig.suptitle(f'N={N},T1={T1},T2={T2}',fontsize = 16)
            self.fig.tight_layout(rect=[0, 0, 1, 0.95])
            self.fig.savefig(f"{chemin_image}")
            self.fig.suptitle('')
            messagebox.showinfo('Sauvegarder',"Image enregstréé")   

        bouton_plus = tk.Button(self, text="+ de sites", 
                                  command=plus)
        bouton_plus.grid(row=1, column=0, pady=10, sticky="nsew")

        bouton_moins = tk.Button(self, text="- de sites", 
                                  command=moins)
        bouton_moins.grid(row=2, column=0, pady=10, sticky="nsew")

        bouton_enregistrer = tk.Button(self, text="Enregistrer", 
                                  command=enregistrer)
        bouton_enregistrer.grid(row=3, column=0, pady=10, sticky="nsew")

        #sliders
        self.T2_slider = tk.Scale(self, from_=0, to=500, resolution=10,
                              orient='horizontal', label='T2', length=200,
                              command=self.update_plot)
        self.T2_slider.set(50)
        self.T2_slider.grid(row=0, column=2)

        self.T1_slider = tk.Scale(self, from_=0, to=500, resolution=10,
                              orient='horizontal', label='T1', length=200,
                              command=self.update_plot)
        self.T1_slider.set(50)
        self.T1_slider.grid(row=0, column=1)

        self.force_index = tk.IntVar()
        self.force_menu = ttk.Combobox(self, textvariable=self.force_index, values=self.force_index, width=5)
        #self.force_menu.grid(row=3, column=0)

        # Figure matplotlib

        self.canvas = None 
        self.fig = None
        self.afficher_subplots()  # Crée self.fig et self.axs
        self.update_plot(None)

        # Configurer les lignes et colonnes pour centrer le contenu
        self.grid_rowconfigure(0, weight=3)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self.grid_rowconfigure(3, weight=1)

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)

    def plot_modes(self,w,amp):
        x=np.arange(1,len(w)+1,)
        if len(w)>12:
            amp=amp[:12,:]
        
        for i in range(min(len(self.axs), amp.shape[0])):
            ax = self.axs[i]
            ax.clear()
            if amp[i, 0] < 0:
                amp[i] *= -1
            ax.bar(x, amp[i], color='skyblue', edgecolor='black',linewidth=1,width=1)
            ax.set_title(f'Mode {i+1}')
            ax.axis('off')
        N = self.centre.Valeurs.get('N', 3)
        if N>11:
            N=11

        ax = self.axs[N]
        ax.clear()
        ax.scatter(x, w, color='skyblue')#, edgecolor='black',)
        ax.set_title('Valeurs propres')
        ax.axis('on')
            
        self.fig.tight_layout()
        self.canvas.draw()

          

    def update_plot(self, event=None):
        if self.update_after_id:
            self.after_cancel(self.update_after_id)

        self.update_after_id = self.after(300, self._do_update_plot)  # 100 ms debounce

    def _do_update_plot(self):
        self.update_after_id = None  # reset debounc
        
        T2 = self.T2_slider.get()
        N = self.centre.Valeurs.get('N', 3)
        l = self.centre.Valeurs.get('longueur', 1)
        T1 = self.T1_slider.get()
        m = self.centre.Valeurs.get('masse', 1)
        self.centre.Valeurs['T2'] = T2
        self.centre.Valeurs['T1'] = T1

        poss=[-T1,-T2] #alternance des données en dessous ou au dessus de la diagonale supérieure
        perte=[g/l + (T1+T2),g/l + (T1+T2)] #alternance des données sur la diagonale principale
        perte=[0,0]
        données=np.zeros([len(poss),2])
        for k in range(len(poss)):
            données[k,0]=poss[k]
            données[k,1]=perte[k]
            limites=(0,0,0)

        w,amp=v_propres(Matrice,données,limites,N)
        indices_trie = np.argsort(w)
        w = w[indices_trie]
        amp = amp.T[indices_trie]
        self.plot_modes(w,amp)

    def afficher_subplots(self):
        N = self.centre.Valeurs.get('N', 2)  # Valeur entrée par l'utilisateur

        if self.canvas:
            self.canvas.get_tk_widget().destroy()  # Supprimer l'ancien widget

        if N>10:
            N=10
        self.fig, axs = plt.subplots(2, ceil(N/2)+1)#, figsize=(min(4 , 12), 4))
        self.axs = axs.flatten() if N > 1 else [axs]
        
        for i, ax in enumerate(self.axs):
            if i < N+2:
                #ax.plot([0, 1], [0, i + 1])
                ax.set_title(" ")
                ax.axis('off')
            else:
                ax.set_visible(False)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(row=1, column=1,rowspan=3, columnspan=3, sticky="nsew")
        self.canvas.draw()

    def tkraise(self, *args, **kwargs):
        super().tkraise(*args, **kwargs)

        N = self.centre.Valeurs.get("N", 5)

        val = list(range(1, N + 1))
        self.force_menu['values'] = val
        self.force_menu.set(val[0] if val else 1)
        self.afficher_subplots()
        self.update_plot(None)

class Comparaison_T2_cercle(tk.Frame):
    def __init__(self, Main, centre):
        super().__init__(Main)
        self.centre = centre

        # Créer le contenu de l'autre page
        label = tk.Label(self, text="Comparaison T2 (pas de force externe)", font=("Helvetica", 20))
        label.grid(row=0, column=1, pady=14, sticky="nsew")

        bouton_retour = tk.Button(self, text="Retour en arrière", 
                                  command=lambda: centre.show_page(Pendule), font=("Helvetica", 14))
        bouton_retour.grid(row=0, column=0, pady=10, sticky="nsew")

        # Fonction pour récupérer les valeurs et fermer la fenêtre
        def valider():
            try:
                self.centre.Valeurs['N'] = int(self.entry_sites.get())
                self.centre.Valeurs['longueur'] = float(self.entry_longueur.get())
                #self.centre.Valeurs['T1'] = float(self.entry_T1.get())
                self.centre.Valeurs['masse'] = float(self.entry_masse.get())
                centre.show_page(Graph_comp_T2_cercle)
            
            except ValueError:
                messagebox.showerror("Erreur", "Veuillez entrer des valeurs numériques valides.")

        # Création des widgets
        self.label_sites = tk.Label(self, text="Nombres de sites :")
        self.entry_sites = tk.Entry(self)
        self.label_longueur = tk.Label(self, text="Longueur (m):")
        self.entry_longueur = tk.Entry(self)
        self.entry_longueur.insert(0,'1.0')
        #self.label_T1 = tk.Label(self, text="T1 (N/m):")
        #self.entry_T1 = tk.Entry(self)
        #self.entry_T1.insert(0,'15.0')
        self.label_masse = tk.Label(self, text="Masse (kg):")
        self.entry_masse = tk.Entry(self)
        self.entry_masse.insert(0,'1.0')
        bouton_valider = tk.Button(self, text="Valider", command=valider)

        # Placement des widgets
        self.label_sites.grid(row=1, column=0, padx=10, pady=5)
        self.entry_sites.grid(row=1, column=1, padx=10, pady=5)
        self.label_longueur.grid(row=2, column=0, padx=10, pady=5)
        self.entry_longueur.grid(row=2, column=1, padx=10, pady=5)
        #self.label_T1.grid(row=3, column=0, padx=10, pady=5)
        #self.entry_T1.grid(row=3, column=1, padx=10, pady=5)
        self.label_masse.grid(row=3, column=0, padx=10, pady=5)
        self.entry_masse.grid(row=3, column=1, padx=10, pady=5)
        bouton_valider.grid(row=4,column=1, pady=10)

        # Configurer les lignes et colonnes pour centrer le contenu
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=2)
        self.grid_rowconfigure(2, weight=2)
        self.grid_rowconfigure(3, weight=2)
        self.grid_rowconfigure(4, weight=2)

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)

class Graph_comp_T2_cercle(tk.Frame):
    def __init__(self, Main, centre):
        super().__init__(Main)
        self.update_after_id = None
        self.centre = centre

        bouton_retour = tk.Button(self, text="Retour en arrière", 
                                  command=lambda: centre.show_page(Comparaison_T2_cercle))
        bouton_retour.grid(row=0, column=0, pady=10, sticky="nsew")

        def plus():
            N = self.centre.Valeurs.get('N', 2)
            N+=2
            self.centre.Valeurs['N']=N
            if N<12:
                self.afficher_subplots()
            self.update_plot(None)
            
        def moins():
            N = self.centre.Valeurs.get('N', 2)
            N-=2
            self.centre.Valeurs['N']=N
            if N<12:
                self.afficher_subplots()
            self.update_plot(None)

        def enregistrer():
            N = self.centre.Valeurs.get('N', 2)
            T1 = self.centre.Valeurs.get('T1', 2)
            T2 = self.centre.Valeurs.get('T2', 2)
            chemin = f'C:/Users/Mathieu/.anaconda/STAGE/modes_normaux_1D_cercle/{N}'
            if os.path.exists(chemin)==False:
                os.makedirs(chemin)
            chemin_image = os.path.join(chemin, f"T1={T1},T2={T2}.png")    
            self.fig.suptitle(f'N={N},T1={T1},T2={T2}',fontsize = 16)
            self.fig.tight_layout(rect=[0, 0, 1, 0.95])
            self.fig.savefig(f"{chemin_image}")
            self.fig.suptitle('')  
            messagebox.showinfo('Sauvegarder',"Image enregstréé")   

        bouton_plus = tk.Button(self, text="+ de sites", 
                                  command=plus)
        bouton_plus.grid(row=1, column=0, pady=10, sticky="nsew")

        bouton_moins = tk.Button(self, text="- de sites", 
                                  command=moins)
        bouton_moins.grid(row=2, column=0, pady=10, sticky="nsew")

        bouton_enregistrer = tk.Button(self, text="Enregistrer", 
                                  command=enregistrer)
        bouton_enregistrer.grid(row=3, column=0, pady=10, sticky="nsew")

        #sliders
        self.T2_slider = tk.Scale(self, from_=0, to=500, resolution=10,
                              orient='horizontal', label='T2', length=200,
                              command=self.update_plot)
        self.T2_slider.set(50)
        self.T2_slider.grid(row=0, column=2)

        self.T1_slider = tk.Scale(self, from_=0, to=500, resolution=10,
                              orient='horizontal', label='T1', length=200,
                              command=self.update_plot)
        self.T1_slider.set(50)
        self.T1_slider.grid(row=0, column=1)

        self.force_index = tk.IntVar()
        self.force_menu = ttk.Combobox(self, textvariable=self.force_index, values=self.force_index, width=5)
        #self.force_menu.grid(row=3, column=0)

        # Figure matplotlib

        self.canvas = None 
        self.fig = None
        self.afficher_subplots()  # Crée self.fig et self.axs
        self.update_plot(None)

        # Configurer les lignes et colonnes pour centrer le contenu
        self.grid_rowconfigure(0, weight=3)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self.grid_rowconfigure(3, weight=1)

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)

    def plot_modes(self,w,amp):
        x=np.arange(1,len(w)+1,)
        if len(w)>12:
            amp=amp[:12,:]
        
        for i in range(min(len(self.axs), amp.shape[0])):
            ax = self.axs[i]
            ax.clear()
            if amp[i, 0] < 0:
                amp[i] *= -1
            ax.bar(x, amp[i], color='skyblue', edgecolor='black',linewidth=1,width=1)
            ax.set_title(f'Mode {i+1}')
            ax.axis('off')
        N = self.centre.Valeurs.get('N', 3)
        if N>11:
            N=11

        ax = self.axs[N]
        ax.clear()
        ax.scatter(x, w, color='skyblue')#, edgecolor='black',)
        ax.set_title('Valeurs propres')
        ax.axis('on')
            
        self.fig.tight_layout()
        self.canvas.draw()

          

    def update_plot(self, event=None):
        if self.update_after_id:
            self.after_cancel(self.update_after_id)

        self.update_after_id = self.after(300, self._do_update_plot)  # 100 ms debounce

    def _do_update_plot(self):
        self.update_after_id = None  # reset debounc
        
        T2 = self.T2_slider.get()
        N = self.centre.Valeurs.get('N', 3)
        l = self.centre.Valeurs.get('longueur', 1)
        T1 = self.T1_slider.get()
        m = self.centre.Valeurs.get('masse', 1)
        self.centre.Valeurs['T2'] = T2
        self.centre.Valeurs['T1'] = T1

        poss=[-T1,-T2] #alternance des données en dessous ou au dessus de la diagonale supérieure
        perte=[g/l + (T1+T2),g/l + (T1+T2)] #alternance des données sur la diagonale principale
        perte=[0,0]
        données=np.zeros([len(poss),2])
        for k in range(len(poss)):
            données[k,0]=poss[k]
            données[k,1]=perte[k]
            limites=(0,0,-T2)

        w,amp=v_propres(Matrice,données,limites,N)
        indices_trie = np.argsort(w)
        w = w[indices_trie]
        amp = amp.T[indices_trie]
        self.plot_modes(w,amp)

    def afficher_subplots(self):
        N = self.centre.Valeurs.get('N', 2)  # Valeur entrée par l'utilisateur

        if self.canvas:
            self.canvas.get_tk_widget().destroy()  # Supprimer l'ancien widget

        if N>10:
            N=10
        self.fig, axs = plt.subplots(2, ceil(N/2)+1, figsize=(min(4 , 12), 4))
        self.axs = axs.flatten() if N > 1 else [axs]
        
        for i, ax in enumerate(self.axs):
            if i < N+2:
                #ax.plot([0, 1], [0, i + 1])
                ax.set_title(" ")
                ax.axis('off')
            else:
                ax.set_visible(False)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(row=1, column=1,rowspan=3, columnspan=3, sticky="nsew")
        self.canvas.draw()

    def tkraise(self, *args, **kwargs):
        super().tkraise(*args, **kwargs)

        N = self.centre.Valeurs.get("N", 5)

        val = list(range(1, N + 1))
        self.force_menu['values'] = val
        self.force_menu.set(val[0] if val else 1)
        self.afficher_subplots()
        self.update_plot(None)

class Comparaison_T2_2D(tk.Frame):
    def __init__(self, Main, centre):
        super().__init__(Main)
        self.centre = centre

        # Créer le contenu de l'autre page
        label = tk.Label(self, text="Comparaison T2 2D (pas de force externe)", font=("Helvetica", 20))
        label.grid(row=0, column=1, pady=14, sticky="nsew")

        bouton_retour = tk.Button(self, text="Retour en arrière", 
                                  command=lambda: centre.show_page(Pendule), font=("Helvetica", 14))
        bouton_retour.grid(row=0, column=0, pady=10, sticky="nsew")

        # Fonction pour récupérer les valeurs et fermer la fenêtre
        def valider():
            try:
                self.centre.Valeurs['N'] = int(self.entry_sites.get())
                self.centre.Valeurs['longueur'] = float(self.entry_longueur.get())
                #self.centre.Valeurs['T1'] = float(self.entry_T1.get())
                self.centre.Valeurs['masse'] = float(self.entry_masse.get())
                centre.show_page(Graph_comp_T2_2D)
            
            except ValueError:
                messagebox.showerror("Erreur", "Veuillez entrer des valeurs numériques valides.")

        # Création des widgets
        self.label_sites = tk.Label(self, text="Nombres de sites :")
        self.entry_sites = tk.Entry(self)
        self.entry_sites.insert(0,'12')
        self.label_longueur = tk.Label(self, text="Longueur (m):")
        self.entry_longueur = tk.Entry(self)
        self.entry_longueur.insert(0,'1.0')
        #self.label_T1 = tk.Label(self, text="T1 (N/m):")
        #self.entry_T1 = tk.Entry(self)
        #self.entry_T1.insert(0,'15.0')
        self.label_masse = tk.Label(self, text="Masse (kg):")
        self.entry_masse = tk.Entry(self)
        self.entry_masse.insert(0,'1.0')
        bouton_valider = tk.Button(self, text="Valider", command=valider)

        # Placement des widgets
        self.label_sites.grid(row=1, column=0, padx=10, pady=5)
        self.entry_sites.grid(row=1, column=1, padx=10, pady=5)
        self.label_longueur.grid(row=2, column=0, padx=10, pady=5)
        self.entry_longueur.grid(row=2, column=1, padx=10, pady=5)
        #self.label_T1.grid(row=3, column=0, padx=10, pady=5)
        #self.entry_T1.grid(row=3, column=1, padx=10, pady=5)
        self.label_masse.grid(row=3, column=0, padx=10, pady=5)
        self.entry_masse.grid(row=3, column=1, padx=10, pady=5)
        bouton_valider.grid(row=4,column=1, pady=10)

        # Configurer les lignes et colonnes pour centrer le contenu
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=2)
        self.grid_rowconfigure(2, weight=2)
        self.grid_rowconfigure(3, weight=2)
        self.grid_rowconfigure(4, weight=2)

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)

class Graph_comp_T2_2D(tk.Frame):
    def __init__(self, Main, centre):
        super().__init__(Main)
        self.update_after_id = None
        self.centre = centre

        bouton_retour = tk.Button(self, text="Retour en arrière", 
                                  command=lambda: centre.show_page(Comparaison_T2_2D))
        bouton_retour.grid(row=0, column=0, pady=10, sticky="nsew")

        def plus():
            N = self.centre.Valeurs.get('N', 2)
            N+=2
            self.centre.Valeurs['N']=N
            if N<12:
                self.afficher_subplots()
            self.update_plot(None)
            
        def moins():
            N = self.centre.Valeurs.get('N', 2)
            N-=2
            self.centre.Valeurs['N']=N
            if N<12:
                self.afficher_subplots()
            self.update_plot(None)

        def enregistrer():
            N = self.centre.Valeurs.get('N', 2)
            T1 = self.centre.Valeurs.get('T1', 2)
            T2 = self.centre.Valeurs.get('T2', 2)
            chemin = f'C:/Users/Mathieu/.anaconda/STAGE/modes_normaux_2D_/{N}'
            if os.path.exists(chemin)==False:
                os.makedirs(chemin)
            chemin_image = os.path.join(chemin, f"T1={T1},T2={T2}.png")    
            self.fig.suptitle(f'N={N},T1={T1},T2={T2}',fontsize = 16)
            self.fig.tight_layout(rect=[0, 0, 1, 0.95])
            self.fig.savefig(f"{chemin_image}")
            self.fig.suptitle('')  
            messagebox.showinfo('Sauvegarder',"Image enregstréé")   

        bouton_plus = tk.Button(self, text="+ de sites", 
                                  command=plus)
        bouton_plus.grid(row=1, column=0, pady=10, sticky="nsew")

        bouton_moins = tk.Button(self, text="- de sites", 
                                  command=moins)
        bouton_moins.grid(row=2, column=0, pady=10, sticky="nsew")

        bouton_enregistrer = tk.Button(self, text="Enregistrer", 
                                  command=enregistrer)
        bouton_enregistrer.grid(row=3, column=0, pady=10, sticky="nsew")

        #sliders
        self.T2_slider = tk.Scale(self, from_=0, to=500, resolution=10,
                              orient='horizontal', label='T2', length=200,
                              command=self.update_plot)
        self.T2_slider.set(50)
        self.T2_slider.grid(row=0, column=2)

        self.T1_slider = tk.Scale(self, from_=0, to=500, resolution=10,
                              orient='horizontal', label='T1', length=200,
                              command=self.update_plot)
        self.T1_slider.set(50)
        self.T1_slider.grid(row=0, column=1)

        self.force_index = tk.IntVar()
        self.force_menu = ttk.Combobox(self, textvariable=self.force_index, values=self.force_index, width=5)
        #self.force_menu.grid(row=3, column=0)

        # Figure matplotlib

        self.canvas = None 
        self.fig = None
        self.afficher_subplots()  # Crée self.fig et self.axs
        self.update_plot(None)

        # Configurer les lignes et colonnes pour centrer le contenu
        self.grid_rowconfigure(0, weight=3)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self.grid_rowconfigure(3, weight=1)

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)

    def plot_modes(self,w,amp):
        x=np.arange(1,len(w)+1,)
        if len(w)>12:
            amp=amp[:12,:]
        
        for i in range(min(len(self.axs), amp.shape[0])):
            for j in w:
                if j!=0:
                    ax = self.axs[i]
                    ax.clear()
                    if amp[i, 0] < 0:
                        amp[i] *= -1
                    ax.bar(x, amp[i], color='skyblue', edgecolor='black',linewidth=1,width=1)
                    ax.set_title(f'Mode {i+1}')
                    ax.axis('off')
        N = self.centre.Valeurs.get('N', 3)
        if N>11:
            N=11

        ax = self.axs[N]
        ax.clear()
        ax.scatter(x, w, color='skyblue')#, edgecolor='black',)
        ax.set_title('Valeurs propres')
        ax.axis('on')
            
        self.fig.tight_layout()
        self.canvas.draw()

          

    def update_plot(self, event=None):
        if self.update_after_id:
            self.after_cancel(self.update_after_id)

        self.update_after_id = self.after(300, self._do_update_plot)  # 100 ms debounce

    def _do_update_plot(self):
        self.update_after_id = None  # reset debounc
        
        T2 = self.T2_slider.get()
        N = self.centre.Valeurs.get('N', 3)
        l = self.centre.Valeurs.get('longueur', 1)
        T1 = self.T1_slider.get()
        m = self.centre.Valeurs.get('masse', 1)
        self.centre.Valeurs['T2'] = T2
        self.centre.Valeurs['T1'] = T1
        if N==6:
            Mat=np.zeros((N,N))
            L = [0,-T1,0,-T1,0,-T2]
            L2 = [0,-T1,0,-T2,0,-T1]
            L3 = [0,-T2,0,-T1,0,-T1]
            for i in range(N):
                if i==0 or i==N-1:
                    Mat[i,:]=L
                elif i<3:
                    Mat[i,:]=L2
                else:
                     Mat[i,:]=L3
                L = L[1:] + [L[0]]
                L2 = L2[1:] + [L2[0]]
                L3 = L3[1:] + [L3[0]]
            w,amp = np.linalg.eig((Mat))
            indices_trie = np.argsort(w)
            w = w[indices_trie]
            amp = (amp[:,indices_trie])[:,:]
            self.plot_modes(w,amp.T)
        if N>1:
            Mat = np.zeros((N, N))
            L=[]
            for i in range(N**2):
                if i%2==0:
                    L.append(0)
                elif (i+1)%4==0:
                    L.append(-T2)
                else:
                    L.append(-T1)
            for i in range(N):
                Mat[i,:] = L[0:N]
                L = L[1:] + [L[0]]
            w,amp = np.linalg.eig((Mat))
            indices_trie = np.argsort(w)
            w = w[indices_trie]
            amp = (amp[:,indices_trie])[:,:]
            self.plot_modes(w,amp.T)

    def afficher_subplots(self):
        N = self.centre.Valeurs.get('N', 2)  # Valeur entrée par l'utilisateur

        if self.canvas:
            self.canvas.get_tk_widget().destroy()  # Supprimer l'ancien widget

        if N>10:
            N=10
        self.fig, axs = plt.subplots(2, ceil(N/2)+1, figsize=(min(4 , 12), 4))
        self.axs = axs.flatten() if N > 1 else [axs]
        
        for i, ax in enumerate(self.axs):
            if i < N+2:
                #ax.plot([0, 1], [0, i + 1])
                ax.set_title(" ")
                ax.axis('off')
            else:
                ax.set_visible(False)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(row=1, column=1,rowspan=3, columnspan=3, sticky="nsew")
        self.canvas.draw()

    def tkraise(self, *args, **kwargs):
        super().tkraise(*args, **kwargs)

        N = self.centre.Valeurs.get("N", 5)

        val = list(range(1, N + 1))
        self.force_menu['values'] = val
        self.force_menu.set(val[0] if val else 1)
        self.afficher_subplots()
        self.update_plot(None)
        
class Comparaison_f_ext(tk.Frame):
    def __init__(self, Main, centre):
        super().__init__(Main)
        self.centre = centre

        # Créer le contenu de l'autre page
        label = tk.Label(self, text="Force externe", font=("Helvetica", 20))
        label.grid(row=0, column=1, pady=14, sticky="nsew")

        bouton_retour = tk.Button(self, text="Retour en arrière", 
                                  command=lambda: centre.show_page(Pendule), font=("Helvetica", 14))
        bouton_retour.grid(row=0, column=0, pady=10, sticky="nsew")

        # Fonction pour récupérer les valeurs et fermer la fenêtre
        def valider():
            try:
                self.centre.Valeurs['N'] = int(self.entry_sites.get())
                self.centre.Valeurs['longueur'] = float(self.entry_longueur.get())
                self.centre.Valeurs['F0'] = float(self.entry_F0.get())

                self.centre.Valeurs['masse'] = float(self.entry_masse.get())
                centre.show_page(Graph_f_ext)
            
            except ValueError:
                messagebox.showerror("Erreur", "Veuillez entrer des valeurs numériques valides.")

        # Création des widgets
        self.label_sites = tk.Label(self, text="Nombres de sites :")
        self.entry_sites = tk.Entry(self)
        self.label_longueur = tk.Label(self, text="Longueur (m):")
        self.entry_longueur = tk.Entry(self)
        self.entry_longueur.insert(0,'1.0')
        self.label_F0 = tk.Label(self, text="Force (N):")
        self.entry_F0 = tk.Entry(self)
        self.entry_F0.insert(0,'50.0')
        self.label_masse = tk.Label(self, text="Masse (kg):")
        self.entry_masse = tk.Entry(self)
        self.entry_masse.insert(0,'1.0')
        bouton_valider = tk.Button(self, text="Valider", command=valider)

        # Placement des widgets
        self.label_sites.grid(row=1, column=0, padx=10, pady=5)
        self.entry_sites.grid(row=1, column=1, padx=10, pady=5)
        self.label_longueur.grid(row=2, column=0, padx=10, pady=5)
        self.entry_longueur.grid(row=2, column=1, padx=10, pady=5)
        self.label_F0.grid(row=3, column=0, padx=10, pady=5)
        self.entry_F0.grid(row=3, column=1, padx=10, pady=5)
        self.label_masse.grid(row=4, column=0, padx=10, pady=5)
        self.entry_masse.grid(row=4, column=1, padx=10, pady=5)
        bouton_valider.grid(row=5, column=1, pady=10)

        # Configurer les lignes et colonnes pour centrer le contenu
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)
        self.grid_rowconfigure(3, weight=1)
        self.grid_rowconfigure(4, weight=1)
        self.grid_rowconfigure(5, weight=4)

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)

class Graph_f_ext(tk.Frame):
    def __init__(self, Main, centre):
        super().__init__(Main)
        self.update_after_id = None
        self.centre = centre

        bouton_retour = tk.Button(self, text="Retour en arrière", 
                                  command=lambda: centre.show_page(Comparaison_f_ext))
        bouton_retour.grid(row=0, column=0, pady=10, sticky="nsew")

        #sliders
        self.T1_slider = tk.Scale(self, from_=0, to=50, resolution=2.5,
                              orient='horizontal', label='T1:', length=200,
                              command=self.update_plot)
        self.T1_slider.set(5)
        self.T1_slider.grid(row=0, column=1)

        self.T2_slider = tk.Scale(self, from_=0, to=50, resolution=2.5,
                              orient='horizontal', label='T2:', length=200,
                              command=self.update_plot)
        self.T2_slider.set(5)
        self.T2_slider.grid(row=0, column=2)

        self.gamma_slider = tk.Scale(self, from_=0.0, to=1, resolution=0.05,
                              orient='vertical', label='Gamma:', length=200,
                              command=self.update_plot)
        self.gamma_slider.set(0)
        self.gamma_slider.grid(row=1, column=0)
        
        tk.Label(self, text="Site excité:").grid(row=0, column=3)
        self.force_index = tk.StringVar()
        self.force_menu = ttk.Combobox(self, textvariable=self.force_index, values=[], width=5)
        self.force_menu.grid(row=0, column=4)
        self.force_menu.bind("<<ComboboxSelected>>", self.update_plot)

        # Configurer les lignes et colonnes pour centrer le contenu
        self.grid_rowconfigure(0, weight=2)
        self.grid_rowconfigure(1, weight=1)

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=5)
        self.grid_columnconfigure(2, weight=5)
        self.grid_columnconfigure(3, weight=0)
        self.grid_columnconfigure(4, weight=10)
        
        self.fig, self.axs = plt.subplots(2, 1, figsize=(4, 12))
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(row=1, column=2, columnspan=5, sticky="nsew")
        self.canvas.draw()
        self.update_plot(None)

    def plot_modes(self, F, wd, amp):
        '''
        indices_trie = np.argsort(w)
        w = np.sqrt(w[indices_trie])
        amp=amp.T
        amp = amp[indices_trie,:]
        gamma = self.gamma_slider.get()
        Aq=Calc_amp_modes_f_ext(amp.T, w, f,wd,gamma)
        Ax=amp.T@Aq
        print(w,amp)
        '''
        for ax in self.axs:
            ax.clear()
        ax=self.axs[0]
        for i in range(len(F)):
            '''for j in range(1, len(wd)):
                threshold = 10*F.max()
                diff = np.abs(amp[j, i] - amp[j - 1, i])
                if np.sign(amp[j, i]) != np.sign(amp[j - 1, i]):
                   amp[j-1, i] = np.nan  # Optional: Break the line to avoid sharp jumps
            '''
            ax.plot(wd, 
                     (np.real((amp[:,i]))), 
                 label=f'Masse {i+1})',linestyle='--')
        if len(F)<8:
            ax.legend()
        ax.set_xlim(left=0)
        ax.set_title('Partie réelle')
        ax.set_ylim(bottom=max(-2.5*F.max(),np.nanmin(amp)), top=min(2.5*F.max(),np.nanmax(amp)))
        
        ax=self.axs[1]
        for i in range(len(F)):
            ax.plot(wd,np.imag(amp[:,i]),label=f'Site {i+1}')
        if len(F)<8:
            ax.legend()
        ax.set_xlim(left=0)
        #ax.set_ylim(bottom=max(-2.5*F.max(),np.nanmin(amp)), top=min(2.5*F.max(),np.nanmax(amp)))
        ax.set_title('Partie complexe')
        

        self.canvas.draw()

    def update_plot(self, event=None):
        if self.update_after_id:
            self.after_cancel(self.update_after_id)

        self.update_after_id = self.after(100, self._do_update_plot)  # 100 ms debounce

    def _do_update_plot(self):
        self.update_after_id = None  # reset debounc
        
        T2 = self.T1_slider.get()
        N = self.centre.Valeurs.get('N', 3)
        l = self.centre.Valeurs.get('longueur', 1)
        T1 = self.T2_slider.get()
        m = self.centre.Valeurs.get('masse', 1)
        F0 = self.centre.Valeurs.get('F0', 25)
        gamma = self.gamma_slider.get()

        force = self.force_index.get()
        if force.strip() != "" and force.isdigit():
            num = int(force) - 1
        else:
             return  
        
        #num = int(self.force_index.get()) - 1
        F = np.zeros(N)
        F[num] = F0

        M = np.eye(N)

        poss=[-T1,-T2] #alternance des données en dessous ou au dessus de la diagonale supérieure
        perte=[g/l + (T1+T2),g/l + (T1+T2)] #alternance des données sur la diagonale principale
        #perte=[0,0]
        données=np.zeros([len(poss),2])
        for k in range(len(poss)):
            données[k,0]=poss[k]
            données[k,1]=perte[k]
            limites=(-T1,-T2,0)

        w,amp=v_propres(Matrice,données,limites,N)
        w_max=np.sqrt(w.max())
        K = Matrice(N,N,données,limites)
        if gamma>1:
            wd = np.linspace(0, 1.3*w_max , 10000)
        else:
            wd = np.linspace(0, 1.3*w_max , 100)
        amp = np.zeros((len(wd), N), dtype=complex)

        for i, omega in enumerate(wd):
            try:
                C = 1j*omega*gamma*M
                #if i == 500:
                # Solve in physical space
                X = np.linalg.solve(K- omega**2 * M-C, F)
                amp[i] = X
                print(K- omega**2 * M-C)
                print(amp[i],X)
            except np.linalg.LinAlgError:
                amp[i] = np.nan
        self.plot_modes(F, wd, amp)

    def afficher_subplots(self):
        N = self.centre.Valeurs.get('N', 2)  # Valeur entrée par l'utilisateur

        if self.canvas:
            self.canvas.get_tk_widget().destroy()  # Supprimer l'ancien widget

        if N>12:
            N=12
        self.fig, axs = plt.subplots(1, 2, figsize=(min(4 , 12), 4))
        self.axs = axs.flatten() if N > 1 else [axs]

        for i, ax in enumerate(self.axs):
            if i < N:
                #ax.plot([0, 1], [0, i + 1])
                ax.set_title(" ")
                ax.axis('off')
            else:
                ax.set_visible(False)
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas.get_tk_widget().grid(row=1, column=1, columnspan=3, sticky="nsew")
        self.canvas.draw()

    def tkraise(self, *args, **kwargs):
        super().tkraise(*args, **kwargs)

        N = self.centre.Valeurs.get("N", 5)

        val = list(range(1, N + 1))
        self.force_menu['values'] = val
        self.force_index.set(str(val[0]))
        self.afficher_subplots()
        self.update_plot(None)


class MainApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.withdraw()  # Cache la fenêtre pendant la config

        default_font = tkFont.nametofont("TkDefaultFont")
        default_font.configure(size=14, family="Helvetica")

        # Met à jour aussi les polices TkText, TkMenu si utilisé :
        tkFont.nametofont("TkTextFont").configure(size=14, family="Helvetica")
        tkFont.nametofont("TkMenuFont").configure(size=14, family="Helvetica")

        self.frames = {}
        self.Valeurs = {}

        self.update_idletasks()

        # Configurer la taille de la fenêtre
        largeur = 1600
        hauteur = 900
        ecran_largeur = self.winfo_screenwidth()
        ecran_hauteur = self.winfo_screenheight()
        x = (ecran_largeur - largeur) // 2
        y = (ecran_hauteur - hauteur) // 2
        self.geometry(f"{largeur}x{hauteur}+{x}+{y}")

        self.deiconify()  # Réaffiche après placement
        
        # Créer les pages
        for F in (Menu, Pendule,Pendule_fext, Comparaison_T2,Graph_comp_T2,Graph_f_ext,Comparaison_f_ext,Comparaison_T2_cercle,Graph_comp_T2_cercle,Comparaison_T2_2D,Graph_comp_T2_2D):#,Graph_bandes):
            page = F(self, self)
            self.frames[F] = page
            page.grid(row=0, column=0, sticky="nsew")

        self.show_page(Menu)

        # Configurer les lignes et colonnes de la fenêtre principale pour centrer les pages
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        self.update()



    def show_page(self, page_class):
        frame = self.frames[page_class]
        frame.tkraise()


# Lancer l'interface
if __name__ == "__main__":
    app = MainApp()
    app.mainloop()



