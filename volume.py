from pylab import *
from no import *
# from sympy.geometry import *

# Volume para o caso de triangulos
class Volume2:
    def __init__(self,p1,p2,p3):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        # self.p1 = No(n1,p1)
        # self.p2 = No(n2,p2)
        # self.p3 = No(n3,p3)
        # self.p1 = (n1,p1)
        # self.p2 = (n2,p2)
        # self.p3 = (n3,p3)

        self.faces = [(self.p1,self.p2), (self.p2,self.p3), (self.p3,self.p1)]
        self.nonorthogonality = 0
        self.skewness = 0
        self.fitness = 0
        # self.P1 = Point(p1[0], p1[1])
        # self.P2 = Point(p2[0], p2[1])
        # self.P3 = Point(p3[0], p3[1])
        # self.T = Triangle(self.P1, self.P2, self.P3)

        self.vizinhos = []
    
        self.nome = 0
        self.ficticio = False
        self.temp = 0
        x,y = self.centroide()
        self.P = array([x,y])
        pass

    def calculaFitness(self):
        def mapa(value, leftMin, leftMax, rightMin, rightMax):
            # Figure out how 'wide' each range is
            leftSpan = leftMax - leftMin
            rightSpan = rightMax - rightMin

            # Convert the left range into a 0-1 range (float)
            valueScaled = float(value - leftMin) / float(leftSpan)

            # Convert the 0-1 range into a value in the right range.
            return rightMin + (valueScaled * rightSpan)

        # Calcula a "skewness" como a média para todas as faces
        self.skewness = self.skewness / (float(len(self.vizinhos)))
        # Calcula a função fitness como o mapeamento de "skewness"
        self.fitness = mapa(self.skewness, 0.3, 0, 0, 10)
        if self.fitness < 0:
            self.fitness = 100
        return self.fitness


    def centroide(self):
        """Calcula o centroide do volume
        
        Returns:
            [x,y] -- [coordenadas do centroide]
        """
        x = (self.p1.valor[0]+self.p2.valor[0]+self.p3.valor[0])/3
        y = (self.p1.valor[1]+self.p2.valor[1]+self.p3.valor[1])/3
        return (x,y)

    def normal(self, face):

        ya = face[0].valor[1]
        yb = face[1].valor[1]
        
        xa = face[0].valor[0]
        xb = face[1].valor[0]

        delta_eta = sqrt((xb-xa)**2 + (yb-ya)**2)

        val1 = (yb-ya)/delta_eta
        val2 = -(xb-xa)/delta_eta

        return array([val1,val2])

    def e_eta(self, face):
        ya = face[0].valor[1]
        yb = face[1].valor[1]
        
        xa = face[0].valor[0]
        xb = face[1].valor[0]

        delta_eta = sqrt((xb-xa)**2 + (yb-ya)**2)
        val1 = (xb-xa)/delta_eta
        val2 = (yb-ya)/delta_eta
        return array([val1,val2])

    def e_qsi(self, A):
        xP, yP = self.centroide()
        xA, yA = A.centroide()

        delta_qsi = sqrt((xA-xP)**2 + (yA-yP)**2)

        try:
            val1 = (xA-xP)/delta_qsi
        except:
            pass
        val2 = (yA-yP)/delta_qsi
        return array([val1,val2])

    def difusao_direta(self, A):
        """Retorna o termo de difusão direta referente a esse volume de controle
        
        Arguments:
            n_face {numero} -- [número correspondente a face]
            A {[volume]} -- [volume de controle vizinho que compartilha a mesma face]
        
        Returns:
            [Di] -- [Difusão Direta]
        """

        for face in self.faces:
            p1 = face[0]
            p2 = face[1]
            for face2 in A.faces:
                p1_A = face2[0]
                p2_A = face2[1]
                if (p1==p1_A and p2==p2_A) or (p2==p1_A and p1==p2_A):
                    face_comum = face

        n = self.normal(face_comum)
        e_qsi = self.e_qsi(A)

        a = face_comum[0].valor
        b = face_comum[1].valor
        # area = self.area()
        ya = a[1]
        yb = b[1]
        xa = a[0]
        xb = b[0]
        area = delta_eta = sqrt((xb-xa)**2 + (yb-ya)**2)

        xP, yP = self.centroide()
        xA, yA = A.centroide()

        delta_qsi = sqrt((xA-xP)**2 + (yA-yP)**2)        


        di = (n.dot(n)/(n.dot(e_qsi))) * (area/delta_qsi)
        return di

    def area(self):
        """Retorna a área desse volume de controle
        
        Returns:
            [area] -- [área do volume de controle]
        """

        xa = self.p1.valor[0]
        xb = self.p2.valor[0]
        xc = self.p3.valor[0]
        ya = self.p1.valor[1]
        yb = self.p2.valor[1]
        yc = self.p3.valor[1]

        M = [[xa,xb,xc],[ya,yb,yc],[1,1,1]]
        area = 0.5*linalg.det(M)
        return area

    def angle_between(self, v1, v2):
        """ Returns the angle in radians between vectors 'v1' and 'v2'::

                >>> angle_between((1, 0, 0), (0, 1, 0))
                1.5707963267948966
                >>> angle_between((1, 0, 0), (1, 0, 0))
                0.0
                >>> angle_between((1, 0, 0), (-1, 0, 0))
                3.141592653589793
        """
        v1_u = v1
        v2_u = v2
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    
    def difusao_cruzada(self, A, volumes):
        """Retorna o termo de difusão cruzada do Volume de Controle
        
        Arguments:
            A {[volume]} -- [volume vizinho]
        
        Returns:
            [número] -- [difusão cruzada]
        """
        for face in self.faces:
            p1 = face[0].nome
            p2 = face[1].nome
            for face2 in A.faces:
                p1_A = face2[0].nome
                p2_A = face2[1].nome
                if (p1==p1_A and p2==p2_A) or (p2==p1_A and p1==p2_A):
                    face_comum = face

        n = self.normal(face_comum)

        e_qsi = self.e_qsi(A)
        e_eta = self.e_eta(face_comum)

        # Aqui calcula a não-ortogonalidade
        # que é o angulo entre n e e_qsi
        aN = self.angle_between(n, e_qsi)
        self.nonorthogonality += aN

        a = face_comum[0].valor
        a = np.delete(a, -1)
        b = face_comum[1].valor
        b = np.delete(b, -1)

        # Aqui calcula a assimetria (skewness)
        # vetor m:
        def perp( a ) :
            b = empty_like(a)
            b[0] = -a[1]
            b[1] = a[0]
            return b
        def seg_intersect(a1,a2, b1,b2) :
            da = a2-a1
            db = b2-b1
            dp = a1-b1
            dap = perp(da)
            denom = dot( dap, db)
            num = dot( dap, dp )
            val = (num / denom.astype(float))*db + b1
            return val
        
        fi = seg_intersect(self.P, A.P, a, b)
        f = (asarray(a) + asarray(b))/2
        m = f-fi
        d = A.P - self.P
        self.skewness += linalg.norm(m)/linalg.norm(d)

        a = np.append(a, 0)
        b = np.append(b, 0)
        Tb = []
        Ta = []
        # Para calcular Tb percorre todos os volumes para ver quem tem ele como um de seus vértices e entao faz a média
        for vol in volumes:
            lista = [vol.p1.valor, vol.p2.valor, vol.p3.valor]
            if array_equal(b,lista[0]) or array_equal(b,lista[1]) or array_equal(b,lista[2]):
                Tb.append(vol.temp)
            if array_equal(a,lista[0]) or array_equal(a,lista[1]) or array_equal(a,lista[2]):
                Ta.append(vol.temp)
            
        Tb = sum(Tb) / len(Tb)
        Ta = sum(Ta) / len(Ta)

        Sdc = -(e_qsi.dot(e_eta))/(n.dot(e_qsi)) * (Tb-Ta)
        return Sdc