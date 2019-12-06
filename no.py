from pylab import *

class No:
    def __init__(self, nome, valor):
        self.nome = nome
        self.valor = asarray(valor)
        self.fronteira = False

        self.vizinhos = []

        self.xl = 0
        self.yl = 0