# import the Network class
from Network_Metric_Analyzer import Network

if __name__ == "__main__":
    # create a NetworkX graph object when the file is given in txt format (eg: DREB1A.txt)
    network = Network.create_Graph()
    print(network)

    # create objects
    DREB1A_obj = Network(network)

    # get the number of proteins and interactions in the PPI network.
    proteins_interactions = DREB1A_obj.num_of_proteins()
    print("Number of proteins and interactions:", proteins_interactions)

    # get the degree of a given protein in the PPI network. (eg: ERF24)
    degree = DREB1A_obj.degree_of_protein()

    # plot the degree distribution as a histogram
    degree_distribution = DREB1A_obj.degree_distribution()
    degree_distribution

    # shortest path between two given proteins of a PPI network.(eg: source:-YAB2 & target:-TIFY8)
    shortest_path = DREB1A_obj.shortest_path()

    # calculate the network diameter of a given PPI network.
    network_diameter = DREB1A_obj.network_diameter()
    print("Diameter of the PPI network:", network_diameter)

    # calculate the clustering coefficient of a given protein in a PPI network. (eg: ERF24)
    clustering_coef = DREB1A_obj.clustering_coef()

    # get the betweenness centrality of a given protein in a PPI network. (eg: ERF24)
    betweenness_cent = DREB1A_obj.betweenness_cent()

    # generate betweenness centrality distribution histogram for a given PPI network.
    betweenness_distribution = DREB1A_obj.betweenness_distribution()
    betweenness_distribution

    # compare two PPI networks
    # create a graph object when the file is given in txt format (eg: DREB2A.txt)
    network2 = Network.create_Graph()
    print(network2)

    # compare two graphs using compare_networks method
    Network.compare_networks(network,network2)




