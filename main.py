from seidel import seidel2d
import pandas as pd
from malha import *


def analitica(x, y):
    return (sin((pi*x)/2)*sin((pi*y)/2))


def main():
    malha = Malha()
    malha.melhora_malha()
    malha.resolve()
    malha.plotar('r')
    temp_analiticas = []
    temp_numericas = []
    diferencas = []
    for vol in malha.volumes:
        x, y = vol.centroide()
        t_analitica = analitica(x, y)
        temp_analiticas.append(t_analitica)
        t_numerica = vol.temp
        temp_numericas.append(t_numerica)
        diff = abs(t_analitica-t_numerica)
        diferencas.append(diff)
    # print("analitica={}, numerica={}, diferenca={}".format(t_analitica,t_numerica,diff))

    diferencas.sort(reverse=True)
    print("Maiores 5 diferencas são: ", diferencas[:5])


if __name__ == '__main__':
    main()
