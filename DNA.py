import random

class DNA:
    def __init__(self, nos, genes=None):
        self.nos = nos
        num_nos = len(nos)
        self.genes = []
        self.fitness = 0
        if genes is not None:
            self.genes = genes
        else:
            for i in range(0,num_nos):
                no = nos[i]
                # n = random.random()
                if no.fronteira==False:
                    num1 = random.uniform(-0.1, 0.01)
                    while num1==0:
                        num1 = random.uniform(-0.1, 0.01)
                    num2 = random.uniform(-0.01, 0.01)
                    while num2==0:
                        num2 = random.uniform(-0.01, 0.01)
                    # num1=num2=0.001
                else:
                    num1=num2=0
                self.genes.append([no.valor[0]+num1,no.valor[1]+num2])
                # self.genes.append([no.valor[0]+random.uniform(0, 0),no.valor[1]+random.uniform(0, 0)])


    def crossover(self, partner):
        newgenes = []
        # num_nos = len(self.genes)
        mid = random.randint(0,len(self.genes))
        # mid = int(len(self.genes)/2)
        for i in range(0, len(self.genes)):
            if i>mid:
                newgenes.append(self.genes[i])
            else:
                newgenes.append(partner.genes[i])
        return DNA(self.nos, newgenes)

    def mutation(self):
        for count in range(0,len(self.genes)):
            if random.random() < 0.1:
                gene = self.genes[count]
                if random.random() < 0.5:
                    if random.random()<0.5:
                        gene[0] += 0.01
                    else:
                        gene[0] -= 0.01
                else:
                    if random.random()<0.5:
                        gene[1] += 0.01
                    else:
                        gene[1] -= 0.01
                # gene[0] += random.uniform(-0.1, 0.1)
                # gene[1] += random.uniform(-0.1, 0.1)
