import triangle as tr
from volume import *
import itertools
from DNA import *
from no import *
import meshio
import optimesh
import random
from numpy import linalg as LA
import math
# from gurobipy import *
from numpy import array
from matplotlib.collections import LineCollection
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay, delaunay_plot_2d

def read_poly(file_name):
    """
    Simple poly-file reader, that creates a python dictionary 
    with information about vertices, edges and holes.
    It assumes that vertices have no attributes or boundary markers.
    It assumes that edges have no boundary markers.
    No regional attributes or area constraints are parsed.
    """

    output = {'vertices': None, 'holes': None, 'segments': None}
    
    # open file and store lines in a list
    file = open(file_name, 'r')
    lines = file.readlines()
    file.close()
    lines = [x.strip('\n').split() for x in lines]
    
    # Store vertices
    vertices= []
    N_vertices, dimension, attr, bdry_markers = [int(x) for x in lines[0]]
    # We assume attr = bdrt_markers = 0
    for k in range(N_vertices):
        label, x, y = [items for items in lines[k+1]]
        vertices.append([float(x), float(y)])
    output['vertices']=array(vertices)
    
    # Store segments
    segments = []
    N_segments, bdry_markers = [int(x) for x in lines[N_vertices+1]]
    for k in range(N_segments):
        label, pointer_1, pointer_2 = [items for items in lines[N_vertices+k+2]]
        segments.append([int(pointer_1)-1, int(pointer_2)-1])
    output['segments'] = array(segments)
    
    # Store holes
    N_holes = int(lines[N_segments+N_vertices+2][0])
    holes = []
    for k in range(N_holes):
        label, x, y = [items for items in lines[N_segments + N_vertices + 3 + k]]
        holes.append([float(x), float(y)])
    
    output['holes'] = array(holes)
    
    return output

class Malha:
    def __init__(self, dna=None):
        self.volumes = []
        self.nos = []
        box1 = tr.get_data('ell')
        '''
        box1 tem a seguinte estrutura:
        triangles -> matriz com todos os triangulos da malha cada triangulo ex: [1,2,3]
        vertex_markers -> marca os vertices que fazem fronteira no dominio
        vertices -> tabela com as coordenadas de todos os vertices da malha
        '''
        self.box2 = tr.triangulate(box1, 'ra0.3')
        '''
        self.box2 tem a mesma estrutura que box1 mas tem a triangulação mais acentuada
        '''

        # la = tr.get_data('la')
        # self.box2 = tr.triangulate(la, 'pq')
        # lake_superior = read_poly("superior.poly")
        # vertices_ls = lake_superior['vertices']
        # tri = Delaunay(vertices_ls)
        # for i in range(len(self.box2['vertices'])):
        #     if self.box2['vertex_markers'][i]==0:
        #         self.box2['vertices'][i][0] += 0.1
        #         self.box2['vertices'][i][1] += 0.1

        self.vols = self.box2['triangles'].tolist()
        self.verts = self.box2['vertices'].tolist()

        # Diminui a escala do domínio para de 0 até 1 em x e y
        for vert in self.verts:
            vert[0] = vert[0]/4
            vert[1] = vert[1]/4

        # Desordena a malha
        self.verts[7][1] -= 0.11
        self.verts[16][0] -= 0.1
        self.verts[16][1] += 0.1
        self.verts[6][0] += 0.1
        self.verts[6][1] += 0.1
        self.verts[8][0] += 0.1
        self.verts[8][1] += 0.1
        self.verts[16][0] += 0.1
        self.verts[16][1] += 0.1

        self.box2['vertices'] = self.verts

        self.vertices_reais = self.verts.copy()
        self.volumes_reais = self.vols.copy()
        self.salva_malha('foo.vtk')
        self.carrega_malha('foo.vtk')
    
    def melhora_malha(self):
        """
        Roda algum algoritmo de melhoramento da malha
        """
        X = asarray(self.vertices_reais)
        cells = asarray(self.volumes_reais)
        # X, cells = optimesh.cpt.fixed_point_uniform(X, cells, Infinity , 50)
        X, cells = optimesh.odt.fixed_point_uniform(X, cells, Infinity, 10000, 1)
        # X, cells = optimesh.cpt.quasi_newton_uniform(X, cells, Infinity, 50)
        self.verts = X
        self.vols = cells

        # self.angular_smoothing()
        # self.verts = X

        # Salva a malha melhorada no arquivo out.vtk
        self.salva_malha('out.vtk')
        self.carrega_malha('out.vtk')
        return

    def carrega_malha(self, file):
        '''
        Baseado na malha foo.vtk, cria os objetos correspondentes
        '''
        self.volumes.clear()
        self.nos.clear()
        # self.verts.clear()
        # self.vols.clear()

        mesh = meshio.read(file)
        self.verts = mesh.points
        self.vols = mesh.cells['triangle']

        cont_vertices = len(self.verts)-1
        
        for count in range(0,len(self.verts)):
            novoNo = No(count, self.verts[count])
            # Define os vizinhos desse novo nó:
            for vol in self.vols:
                if count in vol:
                    for i in range(0,3):
                        if vol[i] != count and vol[i] not in novoNo.vizinhos:
                            novoNo.vizinhos.append(vol[i])

            # if self.verts[count][0] == 0 or self.verts[count][0] == 1 or self.verts[count][1] == 0 or self.verts[count][1] == 1 or (self.verts[count][0]==0.5 and self.verts[count][1]>=0.5) or (self.verts[count][1]==0.5 and self.verts[count][0]>=0.5):
            #     novoNo.fronteira=True
            if self.box2['vertex_markers'][count]==1:
                novoNo.fronteira=True
            
            self.nos.append(novoNo)
        
        for no in self.nos:
            vizinhos = []
            for vizinho_indice in no.vizinhos:
                vizinhos.append(self.nos[vizinho_indice])
            no.vizinhos = vizinhos

        
        # ALGORITMO: Move nós para centróides usando o Gurobi
        # try:
        #     m = Model("modelo2")
        #     obj = 0
        #     for no in self.nos:
        #         # xl = m.addVar(name="xl")
        #         # yl = m.addVar(name="yl")
        #         no.xl = m.addVar(name="xl")
        #         no.yl = m.addVar(name="yl")
        #         # obj = 0
        #         if no.fronteira==False:
        #             Ni = no.valor
        #             # Nj = no.vizinhos[i].valor
        #             # Nj = no.vizinhos[i]
                    
        #             for i in range(0,len(no.vizinhos)):
        #                 # Nj = no.vizinhos[i].valor
        #                 Nj = no.vizinhos[i]
        #                 # Nj_mais1 = no.vizinhos[i+1].valor
        #                 # Nj_menos1 = no.vizinhos[i-1].valor
        #                 # if Nj.xl.X!=0 or Nj.yl.X!=0:
        #                 #     xi = Nj.xl
        #                 #     yi = Nj.yl
        #                 # else:
        #                 xi = Nj.valor[0]
        #                 yi = Nj.valor[1]
        #                 novo_Ni = np.array([no.xl, no.yl])
        #                 obj = obj + (no.xl-xi)*(no.xl-xi)+(no.yl-yi)*(no.yl-yi)
        #     m.setObjective(obj, GRB.MINIMIZE)
        #     m.optimize()
        #     print('obj: %g' % m.objVal)
        #     for no in self.nos:
        #         if no.xl.X!=0 or no.yl.X!=0:
        #             no.valor = np.array([no.xl.X, no.yl.X, 0])
        #             pass
        #     # no.valor = np.array([xl.X,yl.X,0])
        #     pass
        # except GurobiError as e:
        #     print('Error code ' + str(e.errno) + ": " + str(e))

        # except AttributeError:
        #     print('Encountered an attribute error')                

        
        # ALGORITMO: Melhoramento por Angulos ***
        # for no in self.nos:
        #     if no.fronteira==False:
        #         valoresX = []
        #         valoresY = []
        #         for i in range(1,len(no.vizinhos)-1):
        #             Ni = no.valor
        #             Nj_mais1 = no.vizinhos[i+1].valor
        #             Nj_menos1 = no.vizinhos[i-1].valor
        #             Nj = no.vizinhos[i].valor

        #             vj = np.array([Ni[0]-Nj[0],Ni[1]-Nj[1]])
        #             vj_mais1 = np.array([Nj_mais1[0]-Nj[0],Nj_mais1[1]-Nj[1]])
        #             vj_menos1 = np.array([Nj_menos1[0]-Nj[0],Nj_menos1[1]-Nj[1]])

        #             denominador = (LA.norm(vj, 2))*(LA.norm(vj_mais1, 2))
        #             denominador2 = (LA.norm(vj, 2))*(LA.norm(vj_menos1, 2))
        #             # print(denominador)

        #             alpha1 = math.acos(np.dot(vj,vj_mais1)/denominador)
        #             alpha2 = math.acos(np.dot(vj,vj_menos1)/denominador2)

        #             Betaj = (alpha2-alpha1)/2

        #             x, y = Ni[0],Ni[1]
        #             x0, y0 = Nj[0], Nj[1]

        #             xl = x0 + (x-x0)*math.cos(Betaj)-(y-y0)*math.sin(Betaj)
        #             yl = y0 + (x-x0)*math.sin(Betaj)+(y-y0)*math.cos(Betaj)

        #             valoresX.append(xl)
        #             valoresY.append(yl)

        #         xl = sum(valoresX)/len(valoresX)
        #         yl = sum(valoresY)/len(valoresY)

        #         no.valor = np.array([xl,yl,0])


        cont = 0
        for volume in self.vols:
            fronteira = False
            pt1 = volume[0]
            pt2 = volume[1]
            pt3 = volume[2]

            vol = Volume2(self.nos[pt1],self.nos[pt2],self.nos[pt3])

            vol.nome = cont
            cont+=1
            self.volumes.append(vol)

            # Se esse volume fizer fronteira, cria o ficticio correspondente
            
            lista = [pt1,pt2,pt3]
            combinacoes = list(itertools.combinations(lista,2))
            for i in combinacoes:
                ponto1 = i[0]
                ponto2 = i[1]
                if pt1!=ponto1 and pt1!=ponto2:
                    ponto3 = pt1
                if pt2!=ponto1 and pt2!=ponto2:
                    ponto3 = pt2
                if pt3!=ponto1 and pt3!=ponto2:
                    ponto3 = pt3

                v_ponto1 = self.nos[ponto1]
                v_ponto2 = self.nos[ponto2]
                v_ponto3 = self.nos[ponto3]
                # v_ponto1 = self.verts[ponto1]
                # v_ponto2 = self.verts[ponto2]
                # v_ponto3 = self.verts[ponto3]
                # 5) Fronteira topo
                if (v_ponto1.valor[1]==1 and v_ponto2.valor[1]==1):
                    fronteira = True
                    # Cria o nó fictício
                    f_pt1 = v_ponto1
                    f_pt2 = v_ponto2
                    valor = [v_ponto3.valor[0],v_ponto1.valor[1]+(v_ponto1.valor[1]-v_ponto3.valor[1]), 0]
                    cont_vertices+=1
                    f_pt3 = No(cont_vertices, valor)
                    # self.verts.append(valor)
                    self.verts = np.vstack((self.verts, np.asarray(valor)))
                    
                    # Cria o volume fictício
                    vol = Volume2(f_pt1,f_pt2,f_pt3)
                    # vol = Volume2(f_pt1, f_pt2, f_pt3,ponto1,ponto2,cont_vertices)

                    vol.ficticio = True
                    vol.nome = cont
                    cont+=1
                    self.volumes.append(vol)
                    continue
                # 1) Fronteira esquerda
                if (v_ponto1.valor[0]==0 and v_ponto2.valor[0]==0):
                    fronteira = True
                    # Cria volume ficticio na esquerda
                    f_pt1 = v_ponto1
                    f_pt2 = v_ponto2
                    valor = [-v_ponto3.valor[0],v_ponto3.valor[1], 0]
                    cont_vertices+=1
                    f_pt3 = No(cont_vertices, valor)
                    self.verts = np.vstack((self.verts, np.asarray(valor)))
                    # self.verts.append(valor)
                    
                    vol = Volume2(f_pt1, f_pt2, f_pt3)
                    vol.ficticio = True
                    vol.nome = cont

                    cont+=1
                    self.volumes.append(vol)
                    continue
                # 2) Fronteira baixo
                if (v_ponto1.valor[1]==0 and v_ponto2.valor[1]==0):
                    fronteira = True
                    # Cria volume ficticio na esquerda
                    f_pt1 = v_ponto1
                    f_pt2 = v_ponto2
                    valor = [v_ponto3.valor[0],-v_ponto3.valor[1], 0]
                    cont_vertices+=1
                    f_pt3 = No(cont_vertices, valor)
                    # self.verts.append(valor)
                    self.verts = np.vstack((self.verts, np.asarray(valor)))
                    
                    vol = Volume2(f_pt1, f_pt2, f_pt3)
                    vol.ficticio = True
                    vol.nome = cont
                    

                    cont+=1
                    self.volumes.append(vol)
                    continue
                # 3) Fronteira direita
                if (v_ponto1.valor[0]==1 and v_ponto2.valor[0]==1):
                    fronteira = True
                    # Cria volume ficticio na esquerda
                    f_pt1 = v_ponto1
                    f_pt2 = v_ponto2
                    valor = [v_ponto1.valor[0]+(v_ponto1.valor[0]-v_ponto3.valor[0]),v_ponto3.valor[1], 0]
                    cont_vertices+=1
                    f_pt3 = No(cont_vertices, valor)
                    # self.verts.append(valor)
                    self.verts = np.vstack((self.verts, np.asarray(valor)))
                    
                    vol = Volume2(f_pt1, f_pt2, f_pt3)
                    vol.ficticio = True
                    vol.nome = cont
                    cont+=1
                    self.volumes.append(vol)
                    continue
                # 4) Fronteira direita meio
                if (v_ponto1.valor[0]==0.5 and v_ponto1.valor[1]>=0.5 and v_ponto2.valor[0]==0.5 and v_ponto2.valor[1]>=0.5):
                    fronteira = True
                    # Cria volume ficticio na esquerda
                    f_pt1 = v_ponto1
                    f_pt2 = v_ponto2
                    valor = [v_ponto1.valor[0]+(v_ponto1.valor[0]-v_ponto3.valor[0]),v_ponto3.valor[1], 0]
                    
                    if valor[0] > 0.5 and valor[1] > 0.5:
                        cont_vertices+=1
                        f_pt3 = No(cont_vertices, valor)
                        # self.verts.append(valor)
                        self.verts = np.vstack((self.verts, np.asarray(valor)))
                        
                        vol = Volume2(f_pt1, f_pt2, f_pt3)
                        vol.ficticio = True
                        vol.nome = cont
                        cont+=1
                        self.volumes.append(vol)
                    continue
                # 5) Fronteira cima meio
                if (v_ponto1.valor[1]==0.5 and v_ponto1.valor[0]>=0.5 and v_ponto2.valor[1]==0.5 and v_ponto2.valor[0]>=0.5):
                    fronteira = True
                    # Cria volume ficticio na esquerda
                    f_pt1 = v_ponto1
                    f_pt2 = v_ponto2
                    valor = [v_ponto3.valor[0],v_ponto1.valor[1]+(v_ponto1.valor[1]-v_ponto3.valor[1]), 0]
                    
                    if valor[1] > 0.5 and valor[0] > 0.5:
                        cont_vertices+=1
                        f_pt3 = No(cont_vertices, valor)
                        # self.verts.append(valor)
                        self.verts = np.vstack((self.verts, np.asarray(valor)))
                        
                        vol = Volume2(f_pt1, f_pt2, f_pt3)
                        vol.ficticio = True
                        vol.nome = cont
                        cont+=1
                        self.volumes.append(vol)
                    continue

        # Define a vizinhança de cada triangulo
        for vol in self.volumes:
            v1 = vol.p1.nome
            v2 = vol.p2.nome
            v3 = vol.p3.nome
            
            lista = [v1,v2,v3]
            combinacoes = list(itertools.permutations(lista,2))
            for possivel in self.volumes:
                if vol==possivel:
                    continue
                p_v1 = possivel.p1.nome
                p_v2 = possivel.p2.nome
                p_v3 = possivel.p3.nome

                p_lista = [p_v1,p_v2,p_v3]
                p_combinacoes = list(itertools.permutations(p_lista,2))
                for i in combinacoes:
                    for j in p_combinacoes:
                        if i==j:
                            if possivel not in vol.vizinhos:
                                # print("{},{},{} tem vizinho {},{},{}".format(v1,v2,v3,p_v1,p_v2,p_v3))
                                vol.vizinhos.append(possivel)
                            break
        # self.plotar()
        return

    def plotar(self, cor):
        plt.clf()
        # plt.ion()
        # plt.show()
        def analitica(x,y):
            return sin((pi*x)/2)*sin((pi*y)/2)

        for vol in self.volumes:
            x,y=vol.centroide()
            diff = str(round(analitica(x,y),2))+' '+str(round(vol.temp, 2))+':'+str(vol.nome)
            if vol.ficticio==True:
                text = diff
                plt.text(x,y,text)
                plt.plot(x, y, 'bo')
            else:
                plt.plot(x, y, 'go')
                text = diff
                plt.text(x,y,text)
            plt.plot([vol.p1.valor[0],vol.p2.valor[0],vol.p3.valor[0],vol.p1.valor[0]], [vol.p1.valor[1],vol.p2.valor[1],vol.p3.valor[1],vol.p1.valor[1]], cor)
        # plt.figure(figsize=(15, 5))
        # ax1 = plt.subplot(131, aspect='equal')
        # tr.plot(ax1, **self.box2)
        # plt.draw()
        plt.show()
        # plt.pause(0.001)

    def resolve(self, it=10):
        '''
        Resolve o problema especificado na dissertação com a malha atual
        '''
        faces = 0
        N = len(self.volumes)
        a = zeros([N,N])
        b = zeros(N)
        for _ in range(0,it):
            cont = 0
            for vol in self.volumes:
                x,y = vol.centroide()
                Sdc = 0 # Difusao Cruzada
                Dsd = 0 # Difusao Direta
                if vol.nome==15:
                    pass
                if vol.ficticio == True:
                    aP = 1
                    aA = -1
                    vizinho = vol.vizinhos[0] # Volumes fictícios tem apenas 1 vizinho
                    a[cont][vizinho.nome] = -aA

                    # Determina qual o valor para Tm
                    Tm = 0
                    lista = [vol.p1.nome,vol.p2.nome,vol.p3.nome]
                    combinacoes = list(itertools.permutations(lista,2))
                    for i in combinacoes:
                        n_p1 = i[0]
                        n_p2 = i[1]
                        v_p1 = self.verts[n_p1]
                        v_p2 = self.verts[n_p2]
                        if v_p1[1]==1 and v_p2[1]==1:
                            Tm = sin(pi*x/2)
                        if v_p1[0]==0.5 and v_p2[0]==0.5 and v_p1[1]>=0.5 and v_p2[1]>=0.5:
                            Tm = (sqrt(2)/2)*sin(pi*x/2)
                        if v_p1[1]==0.5 and v_p2[1]==0.5 and v_p1[0]>=0.5 and v_p2[0]>=0.5:
                            Tm = sin(pi*y/2)

                    bP = 2*Tm

                else:

                    for vizinho in vol.vizinhos:
                        Sdc += vol.difusao_cruzada(vizinho, self.volumes)
                        faces += 1
                        direta = vol.difusao_direta(vizinho)
                        a[cont][vizinho.nome] = -direta
                        Dsd += direta
                    aP = Dsd
                    
                    # Cálculo de bP
                    Sp = -((pi**2)/2)*sin((pi*x)/2)*sin((pi*y)/2)
                    bP = -Sp * vol.area() + Sdc

                a[cont][cont] = aP
                b[cont] = bP
                cont+=1

            T = solve(a,b)
            # Seta os volumes com as temperaturas correspondentes
            for cont in range(0,N):
                vol = self.volumes[cont]
                vol.temp = T[cont]

        # self.nonorthogonality = self.nonorthogonality/faces
        # self.skewness = self.skewness/faces
        # self.calcula_fitness()
        # return T
        return

    # def calcula_fitness(self):
    #     cont = 0
    #     for vol in self.volumes:
    #         if vol.ficticio==False:
    #             vol.calculaFitness()
    #             self.fitness += abs(vol.fitness)
    #             cont += 1
    #     self.fitness = self.fitness / cont

    def salva_malha(self, file):
        ''' Salva a malha atual no formato VTK'''
        points = np.asarray(self.verts)
        cells = {
            "triangle": np.asarray(self.vols)
        }
        meshio.write_points_cells(
            file,
            points,
            cells,
            )
        return