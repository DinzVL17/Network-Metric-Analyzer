'''
Title: Network metric analyzer to analyze a given protein-protein interaction (PPI) network.
Author: Dinuri Vishara Lokuwalpola
Date: 24/04/2023
'''

# import necessary packages
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

# create a class for PPI networks
class Network:

    # constructor method
    def __init__(self, ppi_network):
        self.ppi_network = ppi_network


    # method to create and return a Networkx graph object when a PPI network file
    # from the STRING database is given as an input in txt format.
    @staticmethod
    def create_Graph():
        # get the PPI network file as a user input
        file = input("Enter PPI network file:")
        with open(file, 'r') as networkFile:
            pn = nx.Graph()
            for lines in networkFile:
                lines = lines.strip()
                if lines != "\n":
                    if "#" != lines[0]:
                        lines = lines.split("\t")
                        pn.add_edge(lines[0], lines[1], weight=lines[-1])
            return (pn)

    # method to return the number of proteins and interactions in a given PPI network
    def num_of_proteins(self):
        nodes = self.ppi_network.number_of_nodes()
        edges = self.ppi_network.number_of_edges()

        return nodes, edges

    # a method to calculate the degree of a given protein in a given PPI network.
    def degree_of_protein(self):
        protein = input("Enter a protein:")
        if protein not in list(self.ppi_network.nodes):
            print("Invalid node")
        else:
            degree = self.ppi_network.degree(protein)
            print(f"Degree of {protein} protein: {degree}")

    # method to plot the degree distribution as a histogram of a given PPI network.
    def degree_distribution(self):
        degreeDist = self.ppi_network.degree()
        sns.histplot(degreeDist, bins=30, kde=True)
        plt.xlabel("Degrees of nodes")
        plt.ylabel("Number of proteins")
        plt.title("Histogram of degree distribution")
        plt.show()

    # method to return the shortest path between two given proteins of a PPI network.
    def shortest_path(self):
        source = input("Enter a source protein:")
        target = input("Enter a target protein:")
        if source not in list(self.ppi_network.nodes):
            print("Invalid source node")
        elif target not in list(self.ppi_network.nodes):
            print("Invalid target node")
        else:
            sPath = nx.shortest_path(self.ppi_network,source,target)
            print(f"Shortest path between {source} and {target}: {sPath}")


    # method to calculate the network diameter of a given PPI network.
    def network_diameter(self):
        diameter = nx.diameter(self.ppi_network)
        return diameter

    # method to return the clustering coefficient of a given protein in a PPI network.
    def clustering_coef(self):
        protein = input("Enter a protein:")
        if protein not in list(self.ppi_network.nodes):
            print("Invalid node")
        else:
            cCoefficient = nx.clustering(self.ppi_network, protein)
            print(f"Clustering coefficient of {protein} protein: {cCoefficient}")

    # method to return the betweenness centrality of a given protein in a PPI network.
    def betweenness_cent(self):
        protein = input("Enter a protein:")
        if protein not in list(self.ppi_network.nodes):
            print("Invalid node")
        else:
            # create a dictionary of betweenness centrality
            bc_diction = nx.betweenness_centrality(self.ppi_network)
            # get the betweenness centrality of the given protein
            for value in bc_diction:
                if value == protein:
                    bc = bc_diction[value]
                    print(f"Betweenness centrality of {protein} protein: {bc}")

    # method to generate a betweenness centrality distribution
    # histogram for a given PPI network.
    def betweenness_distribution(self):
        # create a dictionary of betweenness centrality
        bc_diction = nx.betweenness_centrality(self.ppi_network)
        # create the histogram
        sns.histplot(bc_diction, bins=30, kde=True)
        plt.xlabel("Betweenness centralities of nodes")
        plt.ylabel("Number of proteins")
        plt.title("Betweenness centrality distribution ")
        plt.show()

    # method to compare two given PPI networks
    @staticmethod
    def compare_networks(network1,network2):
        # calculate the average clustering coefficients of two networks
        cCoefficient1 = nx.average_clustering(network1)
        cCoefficient2 = nx.average_clustering(network2)
        print("Average clustering coefficient of network 1:", cCoefficient1)
        print("Average clustering coefficient of network 2:", cCoefficient2)

        # calculate the network diameters for two networks
        diameter1 = nx.diameter(network1)
        diameter2 = nx.diameter(network2)
        print("Diameter of network 1:", diameter1)
        print("Diameter of network 2:", diameter2)

        # compare the values using bar plots
        fig, axs = plt.subplots(1, 2, figsize=(10, 4))
        axs[0].bar(["network1", "network2"], [cCoefficient1, cCoefficient2])
        axs[0].set_title("Average clustering coefficients of two networks")
        axs[1].bar(["network1", "network2"], [diameter1, diameter2])
        axs[1].set_title("Network diameters for two networks")
        plt.show()

        # get the degree distributions of two networks.
        degree1 = dict(network1.degree).values()
        degree2 = dict(network2.degree).values()

        # get the betweenness centrality distributions of two networks.
        bcentrality1 = dict(nx.betweenness_centrality(network1)).values()
        bcentrality2 = dict(nx.betweenness_centrality(network2)).values()

        # compare the degree distributions and betweenness centrality of two networks using histograms.
        fig, axs = plt.subplots(1, 2, figsize=(12, 4))
        axs[0].hist(degree1, label='network1', alpha=0.9, edgecolor='red')
        axs[0].hist(degree2, label='network2', alpha=0.5, edgecolor='yellow')
        axs[0].legend()
        axs[0].set(title="Degree distribution of two networks", xlabel="Degrees of nodes")
        axs[1].hist(bcentrality1, label='network1', alpha=0.9, edgecolor='red')
        axs[1].hist(bcentrality2, label='network2', alpha=0.5, edgecolor='yellow')
        axs[1].legend()
        axs[1].set(title="Betweenness centrality distributions of two networks", xlabel="Betweenness centrality")
        plt.show()


